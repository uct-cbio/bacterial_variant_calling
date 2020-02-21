#!/usr/bin/env nextflow

/*
========================================================================================
                         nf-core/rnaseq
========================================================================================
 nf-core Analysis Pipeline.
 #### Homepage / Documentation
 TBD
----------------------------------------------------------------------------------------


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
        --gff                         Path to GFF3 file
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



Channel.fromPath("$baseDir/assets/where_are_my_files.txt")
       .into{ch_where_trim_galore; ch_where_star; ch_where_hisat2; ch_where_hisat2_sort}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

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


ch_mdsplot_header = Channel.fromPath("$baseDir/assets/mdsplot_header.txt")
ch_heatmap_header = Channel.fromPath("$baseDir/assets/heatmap_header.txt")
ch_biotypes_header = Channel.fromPath("$baseDir/assets/biotypes_header.txt")
Channel.fromPath("$baseDir/assets/where_are_my_files.txt")
       .into{ch_where_trim_galore; ch_where_star; ch_where_hisat2; ch_where_hisat2_sort}

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

// Preset trimming options
if (params.pico){
    clip_r1 = 3
    clip_r2 = 0
    three_prime_clip_r1 = 0
    three_prime_clip_r2 = 3
    forward_stranded = true
    reverse_stranded = false
    unstranded = false
}


//Validate inputs
if ( params.genome == false ) {
    exit 1, "Must set a reference genome fasta file (--genome)"
}

if ( params.reads == false ) {
    exit 1, "Must set the path to the sample file (--reads) in csv format"
}


if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtf_makeSTARindex; gtf_makeHisatSplicesites; gtf_makeHISATindex; gtf_makeBED12;
              gtf_star; gtf_dupradar; gtf_featureCounts; gtf_stringtieFPKM }
} else if( params.gff ){
  gffFile = Channel.fromPath(params.gff)
                   .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
} else {
    exit 1, "No GTF or GFF3 annotation specified!"
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


// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'


/*
 * PREPROCESSING - Convert GFF3 to GTF
 */
if(params.gff){
  process convertGFFtoGTF {
      tag "$gff"

      input:
      file gff from gffFile

      output:
      file "${gff.baseName}.gtf" into gtf_makeSTARindex, gtf_makeBED12,
            gtf_star, gtf_dupradar, gtf_featureCounts

      script:
      """
      gffread $gff -T -o ${gff.baseName}.gtf
      """
  }
}

/*
 * PREPROCESSING - Build BED12 file
 */
if(!params.bed12){
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

        script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}



/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version &> v_ngi_rnaseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    bwa --version &> v_bwa.txt
    preseq &> v_preseq.txt
    read_duplication.py --version &> v_rseqc.txt
    echo \$(bamCoverage --version 2>&1) > v_deeptools.txt
    featureCounts -v &> v_featurecounts.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    samtools --version &> v_samtools.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
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


process '1D_prepare_samples' {

  publishDir "$params.SRAdir", mode: "link"

  input:
      file samples from sample_sheet
  output:
      file "sample_sheet_new.csv" into newSampleSheet
      file "*.fastq" into SRA_new_reads
  script:
  """
  python3 /vcf2fasta/process_samples.py -i $samples -f $params.SRAdir/
  """
}

newSampleSheet
  .splitCsv(header:true)
  .map{ row-> tuple(row.number, file(row.R1), file(row.R2)) }
  .set { newSampleChannel }

newSampleSheet
  .splitCsv(header:true)
  .map{ row-> tuple(row.number, file(row.R1), file(row.R2)) }
  .set { newSampleChannelFastQC }



/*
 * Process 1F: FastQC -  NEED TO EDIT
 */

 process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set number, file(R1), file(R2) from newSampleChannelFastQC

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    cat $R1 $R2 > ${number}_merged.fastq
    fastqc -q ${number}_merged.fastq
    """
}


/*
*process '1E_trim_samples' {
*
*  input:
*    set number, file(R1), file(R2) from newSampleChannel
*  output:
*    file "sample_x_forward_paired.trimmed.fq" into forwardTrimmed
*    file "sample_x_reverse_paired.trimmed.fq" into reverseTrimmed
*    val "$number" into sampleNumber
*  script:
*  """
*  java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 8 -phred33 $R1 $R2 sample_x_forward_paired.trimmed.fq output_forward_unpaired.fq sample_x_reverse_paired.trimmed.fq output_reverse_unpaired.fq SLIDINGWINDOW:4:18 MINLEN:36
*  """
*
*}
*/


/*
 * STEP 2 - Trim Galore! -- TO REPLACE Trimmomatic
 */
process trim_galore {
    label 'low_memory'
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
            else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    set number, file(R1), file(R2) from newSampleChannel
    file wherearemyfiles from ch_where_trim_galore.collect()

    output:
    file "*_R1.fq.gz" into forwardTrimmed
    file "*_R2.fq.gz" into reverseTrimmed
    file "*R1.fq.gz" into forward_trimmed_reads_for_srst2
    file "*R2.fq.gz" into reverse_trimmed_reads_for_srst2
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
    file "where_are_my_files.txt"
    val "$number" into sampleNumber_srst2
    val "$number" into sampleNumber

    script:
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
    if (params.singleEnd) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $R1 $R2
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $R1 $R2
        """
    }
}



/*
 *  END OF PART 1
 *********/



/*
 *
 * Step 1: srst2 (run per sample)  -- edit needed
 * https://github.com/kviljoen/uct-srst2/blob/master/main.nf
 */

process srst2 {
    tag { "srst2.${pairId}" }
    publishDir "${params.outdir}/srst2", mode: "copy"

    input:
    file forward_trimmed_reads_for_srst2
    file reverse_trimmed_reads_for_srst2
    val sampleNumber_srst2

    output:
	file("${sampleNumber_srst2}_srst2*")

    script:
    geneDB = params.gene_db ? "--gene_db $gene_db" : ''
    mlstDB = params.mlst_db ? "--mlst_db $mlst_db" : ''
    mlstdef = params.mlst_db ? "--mlst_definitions $mlst_definitions" : ''
    mlstdelim = params.mlst_db ? "--mlst_delimiter $params.mlst_delimiter" : ''
    """
    srst2 --input_pe $forward_trimmed_reads_for_srst2 $reverse_trimmed_reads_for_srst2 --output ${sampleNumber_srst2}_srst2 --min_coverage $params.min_gene_cov --max_divergence $params.max_gene_divergence $mlstDB $mlstdef $mlstdelim $geneDB
    """
}





/**********
 * PART 2: Mapping
 *
 * Process 2A: Align reads to the reference genome
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
    file "sample_${sampleNumber}_sorted.bai" into bamindexfiles
    file "sample_${sampleNumber}_sorted.bam" into bam_rseqc
    file "sample_${sampleNumber}_sorted.bai" into bamindexfiles_rseqc
    file "sample_${sampleNumber}_sorted.bam" into bam_preseq
    file "sample_${sampleNumber}_sorted.bam" into bam_forSubsamp
    file "sample_${sampleNumber}_sorted.bam" into bam_skipSubsamp
    file "sample_${sampleNumber}_sorted.bam" into bam_featurecounts
    file "*.out" into alignment_logs
  script:
  if( aligner == 'bwa-mem' )
    """
    bwa mem $genome $forwardTrimmed $reverseTrimmed | samtools sort -O BAM -o sample_${sampleNumber}_sorted.bam
    samtools index sample_${sampleNumber}_sorted.bam sample_${sampleNumber}_sorted.bai
    """

  else
    error "Invalid aligner: ${aligner}"

}



/*
 * STEP 4 - RSeQC analysis -- EDIT NEEDED
 */
process rseqc {
    label 'high_memory'
    tag "${bam_rseqc.baseName - '.sorted'}"
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else filename
        }

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam_rseqc
    file index from bamindexfiles_rseqc
    file bed12 from bed_rseqc.collect()

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication
    """
}


/*
 * Step 4.1 Subsample the BAM files if necessary
 */
bam_forSubsamp
    .filter { it.size() > params.subsampFilesizeThreshold }
    .map { [it, params.subsampFilesizeThreshold / it.size() ] }
    .set{ bam_forSubsampFiltered }
bam_skipSubsamp
    .filter { it.size() <= params.subsampFilesizeThreshold }
    .set{ bam_skipSubsampFiltered }

process bam_subsample {
    tag "${bam.baseName - '.sorted'}"

    input:
    set file(bam), val(fraction) from bam_forSubsampFiltered

    output:
    file "*_subsamp.bam" into bam_subsampled

    script:
    """
    samtools view -s $fraction -b $bam | samtools sort -o ${bam.baseName}_subsamp.bam
    """
}


/*
 * Step 4.2 Rseqc genebody_coverage -- edit needed
 */
process genebody_coverage {
    label 'mid_memory'
    tag "${bam.baseName - '.sorted'}"
       publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)           "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)         "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else filename
        }

    when:
    !params.skip_qc && !params.skip_genebody_coverage

    input:
    file bam from bam_subsampled.concat(bam_skipSubsampFiltered)
    file bed12 from bed_genebody_coverage.collect()

    output:
    file "*.{txt,pdf,r}" into genebody_coverage_results

    script:
    """
    samtools index $bam
    geneBody_coverage.py \\
        -i $bam \\
        -o ${bam.baseName}.rseqc \\
        -r $bed12
    mv log.txt ${bam.baseName}.rseqc.log.txt
    """
}





/*
 * STEP 5 - preseq analysis EDIT NEEDED
 */
process preseq {
    tag "${bam_preseq.baseName - '.sorted'}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    when:
    !params.skip_qc && !params.skip_preseq

    input:
    file bam_preseq

    output:
    file "${bam_preseq.baseName}.ccurve.txt" into preseq_results

    script:
    """
    preseq lc_extrap -v -B $bam_preseq -o ${bam_preseq.baseName}.ccurve.txt
    """
}



/*
 * Process 2B: Mark duplicate reads - EDIT for QC
 */

process '2B_mark_duplicates' {
  input:
    file sample_bam from bamfiles
  output:
    file "${sample_bam.baseName}_dedup.bam" into dedup_bamfiles
    file "${sample_bam.baseName}.txt" into dedup_logs
    file "${bam.baseName}.markDups.bam" into bam_md
    file "${bam.baseName}.markDups_metrics.txt" into picard_results
    file "${bam.baseName}.markDups.bam.bai"
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
 * STEP 7 - dupRadar - NEED TO EDIT
 */
process dupradar {
    label 'low_memory'
    tag "${bam_md.baseName - '.sorted.markDups'}"
    publishDir "${params.outdir}/dupradar", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_duprateExpDensCurve.txt") > 0) "scatter_curve_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else "$filename"
        }

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    file bam_md
    file gtf from gtf_dupradar.collect()

    output:
    file "*.{pdf,txt}" into dupradar_results

    script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
    def dupradar_direction = 0
    if (forward_stranded && !unstranded) {
        dupradar_direction = 1
    } else if (reverse_stranded && !unstranded){
        dupradar_direction = 2
    }
    def paired = params.singleEnd ? 'single' :  'paired'
    """
    dupRadar.r $bam_md $gtf $dupradar_direction $paired ${task.cpus}
    """
}



/*
 * STEP 8 Feature counts
 */
process featureCounts {
    label 'low_memory'
    tag "${bam_featurecounts.baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("biotype_counts") > 0) "biotype_counts/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf("_gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    input:
    file bam_featurecounts
    file gtf from gtf_featureCounts.collect()
    file biotypes_header from ch_biotypes_header.collect()

    output:
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt" into geneCounts, featureCounts_to_merge
    file "${bam_featurecounts.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
    file "${bam_featurecounts.baseName}_biotype_counts*mqc.{txt,tsv}" into featureCounts_biotype

    script:
    def featureCounts_direction = 0
    def extraAttributes = params.fcExtraAttributes ? "--extraAttributes ${params.fcExtraAttributes}" : ''
    if (forward_stranded && !unstranded) {
        featureCounts_direction = 1
    } else if (reverse_stranded && !unstranded){
        featureCounts_direction = 2
    }
    // Try to get real sample name
    sample_name = bam_featurecounts.baseName - 'Aligned.sortedByCoord.out'
    """
    featureCounts -a $gtf -g ${params.fcGroupFeatures} -o ${bam_featurecounts.baseName}_gene.featureCounts.txt $extraAttributes -p -s $featureCounts_direction $bam_featurecounts
    featureCounts -a $gtf -g ${params.fcGroupFeaturesType} -o ${bam_featurecounts.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
    cut -f 1,7 ${bam_featurecounts.baseName}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${bam_featurecounts.baseName}_biotype_counts_mqc.txt
    mqc_features_stat.py ${bam_featurecounts.baseName}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${bam_featurecounts.baseName}_biotype_counts_gs_mqc.tsv
    """
}

/*
 * STEP 9 - Merge featurecounts
 */
process merge_featureCounts {
    tag "${input_files[0].baseName - '.sorted'}"
    publishDir "${params.outdir}/featureCounts", mode: 'copy'

    input:
    file input_files from featureCounts_to_merge.collect()

    output:
    file 'merged_gene_counts.txt'

    script:
    //if we only have 1 file, just use cat and pipe output to csvtk. Else join all files first, and then remove unwanted column names.
    def single = input_files instanceof Path ? 1 : input_files.size()
    def merge = (single == 1) ? 'cat' : 'csvtk join -t -f "Geneid,Start,Length,End,Chr,Strand,gene_name"'
    """
    $merge $input_files | csvtk cut -t -f "-Start,-Chr,-End,-Length,-Strand" | sed 's/Aligned.sortedByCoord.out.markDups.bam//g' > merged_gene_counts.txt
    """
}


/*
 *  END OF PART 2
 *********/


/*
 *  Virulence DB stuff
 *********/


if( params.amr_db ) {
	/*
	 * Build resistance database index with Bowtie2
	 */
	process BuildAMRIndex {
		tag { "${amr_db.baseName}" }

		input:
        	file amr_db

        	output:
        	file 'amr.index*' into amr_index

        	"""
        	bowtie2-build $amr_db amr.index --threads ${threads}
		"""
	}

	/*
         * Align reads to resistance database with Bowtie2
         */
	process AMRAlignment {
        	publishDir "${params.out_dir}/Alignment", pattern: "*.bam"

        	tag { dataset_id }

        	input:
        	set dataset_id, file(forward), file(reverse) from amr_read_pairs
        	file index from amr_index.first()

        	output:
        	set dataset_id, file("${dataset_id}_amr_alignment.sam") into amr_sam_files
        	set dataset_id, file("${dataset_id}_amr_alignment.bam") into amr_bam_files

        	"""
        	bowtie2 -p ${threads} -x amr.index -1 $forward -2 $reverse -S ${dataset_id}_amr_alignment.sam
        	samtools view -bS ${dataset_id}_amr_alignment.sam | samtools sort -@ ${threads} -o ${dataset_id}_amr_alignment.bam
        	"""
	}

	process AMRResistome {
        	publishDir "${params.out_dir}/Resistome"

        	tag { dataset_id }

        	input:
        	file amr_db
        	set dataset_id, file(amr_sam) from amr_sam_files

        	output:
        	set dataset_id, file("${dataset_id}_amr_gene_resistome.tsv") into amr_gene_level

        	"""
		csa -ref_fp ${vf_db} -sam_fp ${vf_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_amr_gene_resistome.tsv"
        	"""
	}
}

if( params.vf_db ) {
	/*
         * Build resistance database index with Bowtie2
         */
	process BuildVFIndex {
		tag { "${vf_db.baseName}" }

		input:
        	file vf_db

        	output:
        	file 'vf.index*' into vf_index

        	"""
        	bowtie2-build $vf_db vf.index --threads ${threads}
		"""
	}
	/*
         * Align reads to virulence factor database with Bowtie2
         */
	process VFAlignment {
        	publishDir "${params.out_dir}/Alignment", pattern: "*.bam"

        	tag { dataset_id }

        	input:
        	set dataset_id, file(forward), file(reverse) from vf_read_pairs
        	file index from vf_index.first()

        	output:
        	set dataset_id, file("${dataset_id}_vf_alignment.sam") into vf_sam_files
        	set dataset_id, file("${dataset_id}_vf_alignment.bam") into vf_bam_files

        	"""
        	bowtie2 -p ${threads} -x vf.index -1 $forward -2 $reverse -S ${dataset_id}_vf_alignment.sam
        	samtools view -bS ${dataset_id}_vf_alignment.sam | samtools sort -@ ${threads} -o ${dataset_id}_vf_alignment.bam
        	"""
	}

	process VFResistome {
        	publishDir "${params.out_dir}/Resistome"

        	tag { dataset_id }

        	input:
        	file vf_db
        	set dataset_id, file(vf_sam) from vf_sam_files

        	output:
        	set dataset_id, file("${dataset_id}_vf_gene_resistome.tsv") into vf_gene_level

        	"""
        	csa -ref_fp ${vf_db} -sam_fp ${vf_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_vf_gene_resistome.tsv"
        	"""
	}
}

if( params.plasmid_db ) {
	/*
         * Build plasmid index with Bowtie2
         */
	process BuildPlasmidIndex {
		tag { "${plasmid_db.baseName}" }

		input:
        	file plasmid_db

        	output:
        	file 'plasmid.index*' into plasmid_index

        	"""
        	bowtie2-build $plasmid_db plasmid.index --threads ${threads}
		"""
	}
	/*
         * Align reads to plasmid database with Bowtie2
         */
	process PlasmidAlignment {
        	publishDir "${params.out_dir}/Alignment", pattern: "*.bam"

        	tag { dataset_id }

        	input:
        	set dataset_id, file(forward), file(reverse) from plasmid_read_pairs
        	file index from plasmid_index.first()

        	output:
        	set dataset_id, file("${dataset_id}_plasmid_alignment.sam") into plasmid_sam_files
        	set dataset_id, file("${dataset_id}_plasmid_alignment.bam") into plasmid_bam_files

        	"""
        	bowtie2 -p ${threads} -x plasmid.index -1 $forward -2 $reverse -S ${dataset_id}_plasmid_alignment.sam
        	samtools view -bS ${dataset_id}_plasmid_alignment.sam | samtools sort -@ ${threads} -o ${dataset_id}_plasmid_alignment.bam
        	"""
	}

	process PlasmidResistome {
        	publishDir "${params.out_dir}/Resistome"

        	tag { dataset_id }

        	input:
        	file plasmid_db
        	set dataset_id, file(plasmid_sam) from plasmid_sam_files

        	output:
        	set dataset_id, file("${dataset_id}_plasmid_gene_resistome.tsv") into plasmid_gene_level

        	"""
        	csa -ref_fp ${plasmid_db} -sam_fp ${plasmid_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_plasmid_gene_resistome.tsv"
        	"""
	}
}



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
    file "${f_vcf.baseName}_filtered_snps.recode.vcf" into snp_vcfs_bgzip
  script:
  """
  vcftools --vcf $f_vcf --keep-only-indels --recode --recode-INFO-all --out ${f_vcf.baseName}_filtered_indels
  vcftools --vcf $f_vcf --remove-indels --recode --recode-INFO-all --out ${f_vcf.baseName}_filtered_snps
  """

}

/*
 * Integrate SNPs into reference genome with BCFtools
 */
process BuildConesnsusSequence {
	tag { dataset_id }

	publishDir "${params.out_dir}/Consensus"

	input:
	file snp_vcf_file from snp_vcfs
	file genome from genome_file

	output:
	file("${dataset_id}_consensus.fa") into consensus_files
	file("${dataset_id}_in_list.txt") into ksnp3_configuration

	"""
	bgzip -c $snp_vcf_file
	tabix ${dataset_id}_genome_variants.vcf.gz
	cat $genome | bcftools consensus ${dataset_id}_genome_variants.vcf.gz > ${dataset_id}_consensus.fa
	echo -e "$params.out_dir/Consensus/${dataset_id}_consensus.fa\t$dataset_id" >> ${dataset_id}_in_list.txt
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
    file vcf from snp_vcfs_bgzip
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


if( params.draft ) {
	/*
	 * Create configuration file for kSNP3 using the draft assemblies and user-input reference genome
	 */
	process kSNPDraftAndGenomeConfiguration {
		echo true

		input:
                file draft from draft_genomes

                output:
                file("genome_paths.txt") into genome_config

                shell:
                '''
                #!/bin/sh
                echo "!{genome}\t!{genome.baseName}" > genome_paths.txt
                for d in !{draft};
                do
                        echo "!{draft_path}/${d}\t${d%.*}" >> genome_paths.txt
                done
                '''
	}
}

else {
	/*
	 * Create configuration file for kSNP3 using the user-input reference genome
	 */
	process kSNPGenomeConfiguration {
		echo true

		storeDir 'temporary_files'

		input:
		file genome from genome_file

		output:
		file("genome_paths.txt") into genome_config
		file("$genome") into kchooser_genome

		shell:
		'''
		#!/bin/sh
		base=`echo !{genome} | cut -f1 -d '.'`
		fp=`readlink !{genome}`
		echo "${fp}\t${base}" > genome_paths.txt
		'''
	}
}


/*
 *  END OF PART 4
 *********/



/*
 * STEP 12 MultiQC - EDIT NEEDED
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file multiqc_config from ch_multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('alignment/*') from alignment_logs.collect()
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('rseqc/*') from genebody_coverage_results.collect().ifEmpty([])
    file ('preseq/*') from preseq_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    file ('featureCounts/*') from featureCounts_logs.collect()
    file ('featureCounts_biotype/*') from featureCounts_biotype.collect()
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m picard -m preseq -m rseqc -m featureCounts -m hisat2 -m star -m cutadapt -m fastqc
    """
}