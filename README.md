# ![uct-cbio/bacterial_variant_calling](/assets/cbio_logo.png)

# Bacterial variant calling and phylogenetics
A pipeline for variant calling on bacterial genomes created with Nextflow and singularity / docker

This pipeline is currently under development. If you wish to use it in future, please feel free to watch the repository.

## Quickstart 

    nextflow run uct-cbio/bacterial_variant_calling --reads sample_sheet.csv --genome H37Rv.fa -with-docker bacterial_env

This is assuming you have a sample sheet formatted as described bellow, and a docker image created with VarDock called 'bacterial_env'.

## Basic usage: 
    The typical command for running the pipeline is as follows:
    nextflow run uct-cbio/bacterial_variant_calling --reads sample_sheet.csv --genome refgenome.fa -profile uct_hex
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Path te reference genome against which the reads will be aligned (in fasta format).
      -profile                      Hardware config to use. Currently profile available for UCT's HPC 'uct_hex' - create your own if necessary
    
    Other arguments:
      --outdir                      The output directory where the results will be saved
      --SRAdir                      The directory where reads downloaded from the SRA will be stored
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                     
     
     Example run:
     To run on UCT hex
     1) Start a 'screen' session from the headnode
     2) Start an interactive job using: qsub -I -q UCTlong -l nodes=1:series600:ppn=1 -d `pwd`
     3) A typical command would look something like:

## Sample file
To allow for both local reads and reads from the [SRA](https://www.ncbi.nlm.nih.gov/sra) to be used, the pipeline has the 
ability to pull reads from the SRA based on the accession number (eg, SRR5989977). 

The 'number' column must contain a unique value. 

number | origin | replicate | isolate | R1 | R2
------------ | ------------- | ------------- | ------------- | ------------- | -------------
1 | genomic | 1 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
2 | genomic | 2 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
3 | genomic | 3 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
4 | genomic | 1 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
5 | genomic | 2 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
6 | genomic | 3 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
7 | genomic | 1 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
8 | genomic | 2 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
9 | genomic | 3 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
10 | genomic | 1 | H37Rv | SRR5989977 | 


In the above example, samples 1-9 are locally stored where sample 10 is a control sample from the SRA. 
Including the accession number in the R1 column will result in the reads from the SRA to be downloaded and used in the analysis. 
This must be exported to a csv file, with a comma ',' separating the columns:

    number,origin,replicate,isolate,R1,R2
    1,genomic,1,wgs_sample_1,path/to/reads/reads_R1.fq,path/to/reads/reads_R2.fq
    2,genomic,2,wgs_sample_1,path/to/reads/reads_R1.fq,path/to/reads/reads_R2.fq
    ...
    10,genomic,1,H37Rv,SRR5989977

## Prerequisites
[Nextflow](https://www.nextflow.io), [Docker](https://www.docker.com). All other dependencies are found in the included Docker recipe (VarDock). 

Note: if you are working on UCT hex you can simply use the singularity image specified in the uct_hex profile.

## Documentation

Variant caller used: [Freebayes](https://arxiv.org/abs/1207.3907)

Phylogenetic analysis: [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html)

## Other useful info

[Docker to singularity](https://github.com/singularityware/docker2singularity)


## Built With


## Credits


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
