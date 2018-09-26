

def helpMessage() {
    log.info"""
    =========================================

    Bacterial variant calling and phylogenetics pipeline

    Developed by the bioinformatics support team at the University of Cape Town

    =========================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow main.nf --reads sample_sheet.csv --genome <path to fasta> -with-docker <docker image>

    or

    nextflow main.nf --reads sample_sheet.csv --genome <path to fasta> -with-singularity <singularity image>


    Mandatory arguments:
        --reads                       The sample sheet containing the paths to the fastq files, as well as sample names.
        --genome                      The reference genome to be used in fasta format. Also acts as an outgroup.
        -profile                      Hardware config to use. local / uct_hex

    Optional arguments:
        --minQuality                  The minimum quality to be passed to vcf-tools for filtering variants.
        --vcf_qual_cutoff             Soon to be removed
        --aligner                     Currently only bwa-mem
        --variant_caller              Currently only freebayes


    Other arguments:
        --SRAdir                      The directory where reads downloaded from the SRA will be stored
        --outdir                      The output directory where the results will be saved
        --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.




    """.stripIndent()
}


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Configurable variables
params.name = false
params.project = false
params.email = false
params.plaintext_email = false

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Define the default parameters
 */

//params.outdir            = "$baseDir"

//params.SRAdir            = "$baseDir/"



//Validate inputs
if ( params.genome == false ) {
    exit 1, "Must set a reference genome fasta file (--genome)"
}

if ( params.reads == false ) {
    exit 1, "Must set the path to the sample file (--reads) in csv format"
}


// Has the run name been specified by the user?
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}




log.info """\
Bacterial WGS variant pipeline v0.1
================================
genome   : $params.genome
reads    : $params.reads
Output   : $params.outdir
SRA dir  : $params.SRAdir
"""

/*
 *  Parse the input parameters
 */

genome_file     = file(params.genome)
sample_sheet    = file(params.reads)
reads_ch        = Channel.fromFilePairs(params.reads)
threads         = 4
aligner         = params.aligner
variant_caller  = params.variant_caller
vcf_qual_cutoff = params.vcf_qual_cutoff


/**********
 * PART 1: Data preparation
 *
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */

process '1A_prepare_genome_samtools' { 
  tag "$genome.baseName"
  
  input: 
      file genome from genome_file 
 
  output: 
      file "${genome}.fai" into genome_index_ch  
  
  script:
  """
  samtools faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  PICARD=`which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}

/*
 * Process 1C: Create a FASTA genome sequence dictionary for BWA
 */

process '1C_prepare_genome_bwa' {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "${genome}.amb" into genome_bwa_amb
      file "${genome}.ann" into genome_bwa_ann
      file "${genome}.bwt" into genome_bwa_bwt
      file "${genome}.pac" into genome_bwa_pac
      file "${genome}.sa" into genome_bwa_sa

  script:
  """
  bwa index $genome 
  """
}


process '1D_prepare_samples' {

  publishDir "$params.SRAdir", mode: "link"

  input:
      file samples from sample_sheet
  output:
      file "sample_sheet_new.csv" into newSampleSheet
      file "*.fastq" into SRA_new_reads
  script:
  """
  python3 /vcf2fasta/process_samples.py -i $samples -f $params.SRAdir/ --quick_test_mode True
  """
}

newSampleSheet
  .splitCsv(header:true)
  .map{ row-> tuple(row.number, file(row.R1), file(row.R2)) }
  .set { newSampleChannel }

process '1E_trim_samples' {

  input:
    set number, file(R1), file(R2) from newSampleChannel
  output:
    file "sample_x_forward_paired.trimmed.fq" into forwardTrimmed
    file "sample_x_reverse_paired.trimmed.fq" into reverseTrimmed
    val "$number" into sampleNumber
  script:
  """
  head $R1 > this.txt
  java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 8 -phred33 $R1 $R2 sample_x_forward_paired.trimmed.fq output_forward_unpaired.fq sample_x_reverse_paired.trimmed.fq output_reverse_unpaired.fq SLIDINGWINDOW:4:18 MINLEN:36
  """

}




/*
 *  END OF PART 1
 *********/



/**********
 * PART 2: Mapping
 *
 * Process 2A: Align reads to the genome
 */


/*
 * Process 2A: Map the reads to the reference genome
 */

process '2A_read_mapping' {
  input:
    file forwardTrimmed
    file reverseTrimmed
    val sampleNumber
    file genome from genome_file
    file genome_bwa_amb
    file genome_bwa_ann
    file genome_bwa_bwt
    file genome_bwa_pac
    file genome_bwa_sa
  output:
    file "sample_${sampleNumber}_sorted.bam" into bamfiles
  script:
  if( aligner == 'bwa-mem' )
    """
    bwa mem $genome $forwardTrimmed $reverseTrimmed | samtools sort -O BAM -o sample_${sampleNumber}_sorted.bam
    """
  
  else
    error "Invalid aligner: ${aligner}"
    
}

/*
 * Process 2B: Mark duplicate reads
 */

process '2B_mark_duplicates' {
  input:
    file sample_bam from bamfiles
  output:
    file "${sample_bam.baseName}_dedup.bam" into dedup_bamfiles
    file "${sample_bam.baseName}.txt" into dedup_logs
  script:
    """
    PICARD=`which picard.jar` 
    java -jar \$PICARD MarkDuplicates \
    INPUT=$sample_bam \
    OUTPUT=${sample_bam.baseName}_dedup.bam \
    METRICS_FILE=${sample_bam.baseName}.txt
    """
}


/*
 *  END OF PART 2
 *********/



/**********
 * PART 3: Variant calling
 *
 * Process 3A: Call the variants
 */


process '3A_call_variants' {
  input:
    file genome from genome_file
    file sample_bam from dedup_bamfiles
  output:
    set file("${sample_bam.baseName}.vcf"), file("$sample_bam") into vcf_bam_files


  script:
  if( variant_caller == 'freebayes' )
    """
    freebayes -f $genome -p 1 $sample_bam > need_rename.vcf
    echo "unknown ${sample_bam.baseName}\n" > sample_names.txt
    /bcftools/bcftools reheader need_rename.vcf --samples sample_names.txt -o ${sample_bam.baseName}.vcf

    """
  else if( variant_caller == 'samtools' )
    """
    freebayes -f $genome -p 1 $sample_bam > ${sample_bam.baseName}.vcf
    """
  else
    error "Invalid variant caller: ${variant_caller}"


}

process '3B_calc_coverage' {
  input:
     set file(vcf), file(bam) from vcf_bam_files
  output:
     set file(vcf), file(bam), stdout into vcf_bam_cov_files
  script:
  """
  baseCov=\$(samtools depth $bam | awk '{sum+=\$3} END { print sum/NR}')
  maxCov=\$(echo "\$baseCov * 5" | bc)
  echo \$maxCov
  """
}



process '3C_filter_variants' {
  input:
    set file(vcf), file(bam), coverage from vcf_bam_cov_files
  output:
    file "${vcf.baseName}_filtered.recode.vcf" into filtered_vcfs
  script:
  """
  vcftools --vcf $vcf --minQ $params.minQuality --recode --recode-INFO-all --out ${vcf.baseName}_filtered --maxDP $coverage
  """
}



process '3D_split_vcf_indel_snps' {
  publishDir "${params.outdir}/variants", mode: "link"

  input:
    file f_vcf from filtered_vcfs
  output:
    file "${f_vcf.baseName}_filtered_indels.recode.vcf" into indel_vcfs
    file "${f_vcf.baseName}_filtered_snps.recode.vcf" into snp_vcfs
  script:
  """
  vcftools --vcf $f_vcf --keep-only-indels --recode --recode-INFO-all --out ${f_vcf.baseName}_filtered_indels
  vcftools --vcf $f_vcf --remove-indels --recode --recode-INFO-all --out ${f_vcf.baseName}_filtered_snps
  """

}

/*
 *  END OF PART 3
 *********/



 /**********
 * PART 4: Phylogenetics
 *
 * Process 4A: Create phylogenetic tree
 */


process '4A_bgzip_vcf' {
  input:
    file vcf from snp_vcfs
  output:
    file "${vcf}.gz" into gz_vcfs
    file "${vcf}.gz.tbi" into tbi_vcfs
  script:
  """
  bgzip -c $vcf > ${vcf}.gz
  tabix ${vcf}.gz
  """
}



process '4B_merge_vcf_files' {
  input:
    file gz_vcf from gz_vcfs.collect()
    file tbi_index from tbi_vcfs.collect()
  output:
    file "out_merged.vcf" into merged_vcf
  script:
  """
  /bcftools/bcftools merge $gz_vcf -o out_merged.vcf
  """
}


process '4C_convert_to_phylip_format' {
  input:
    file merged_vcf_file from merged_vcf
  output:
    file "*.phy" into phylip_file
  script:
  """
  python3 /vcf2fasta/vcf2fasta.py -i $merged_vcf_file -o converted -q $vcf_qual_cutoff
  convbioseq -i fasta phylip converted.fa
  """
}


process '4D_run_RAxML' {

  publishDir "${params.outdir}/RAxML", mode: "link", overwrite: false

  input:
    file inphy from phylip_file
    val threads from threads
  output:
    file "*.outFile" into RAxML_out
  
  script:
  """
  /standard-RAxML/raxmlHPC-PTHREADS-AVX -s $inphy -n outFile -m GTRCATX -T $threads -f a -x 123 -N autoMRE -p 456
  """
}

/*
 *  END OF PART 4
 *********/
