

/*
 * Define the default parameters
 */ 

params.genome     = "$baseDir/data/genome.fa"
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.results    = "results"


log.info """\
C A L L I N G S  -  N F    v 1.0 
================================
genome   : $params.genome
reads    : $params.reads
results  : $params.results
"""

/*
 *  Parse the input parameters
 */

genome_file     = file(params.genome)
reads_ch        = Channel.fromFilePairs(params.reads)

designFile = file(params.reads)
reads_ch = Channel
                   .from(designFile)
                   .splitCsv(header: true)
                   .view { sample ->
       println "${sample.number} - ${sample.isolate} - ${sample.R1}"
    }




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

/*
 *  END OF PART 1
 *********/



/**********
 * PART 2: Mapping
 *
 * Process 2A: Align reads to the genome
 */

mode = 'bwa-mem'

process '2_read_mapping' {
  input:
    val sample from reads_ch
    file genome from genome_file
    file genome_bwa_amb
    file genome_bwa_ann
    file genome_bwa_bwt
    file genome_bwa_pac
    file genome_bwa_sa
  output:
    file "sample_${sample.number}_sorted.bam" into bamfiles  
  script:
  if( mode == 'bwa-mem' )
    """
    bwa mem $genome $sample.R1 $sample.R2 | samtools sort -O BAM -o sample_${sample.number}_sorted.bam
    """
  
  else
    error "Invalid alignment mode: ${mode}"
    
}

process 'mark_duplicates' {
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

process 'call_variants' {
  input:
    file genome from genome_file
    file sample_bam from dedup_bamfiles
  output:
    set file("${sample_bam.baseName}.vcf"), file("$sample_bam") into vcf_bam_files
  script:
    """
    freebayes -f $genome -p 1 $sample_bam > ${sample_bam.baseName}.vcf
    """

}

process 'calc_coverage' {
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



process 'filter_variants' {
  input:
    set file(vcf), file(bam), coverage from vcf_bam_cov_files
  output:
    file "${vcf.baseName}_filtered.recode.vcf" into filtered_vcfs
  script:
  """
  vcftools --vcf $vcf --minGQ 20 --recode --recode-INFO-all --out ${vcf.baseName}_filtered --maxDP $coverage
  """
}





