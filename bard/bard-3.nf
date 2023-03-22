#!/usr/bin/env nextflow.

# This pipeline uses Salmon to quantify RNA-Seq data.

# The input data is a FASTQ file containing the reads and a FASTA file containing the reference genome.

# The output data is a BAM file containing the mapped reads and a JSON file containing the quantification results.

# The pipeline is run using the following command:

# nextflow run rnaseq-nf -with-docker

# The pipeline is installed using the following command:

# curl -fsSL get.nextflow.io | bash

# The pipeline is implemented using the following modules:

# - nf-core/rnaseq: This module provides the basic infrastructure for the pipeline.

# - nf-core/salmon: This module provides the Salmon implementation.

# - nf-core/fastqc: This module provides the FastQC implementation.

# - nf-core/umi-tools: This module provides the UMI tools implementation.

# - nf-core/trimgalore: This module provides the Trim Galore! implementation.

# - nf-core/bbsplit: This module provides the BBSplit implementation.

# - nf-core/sortmeRNA: This module provides the SortMeRNA implementation.

# - nf-core/star: This module provides the STAR implementation.

# - nf-core/rsem: This module provides the RSEM implementation.

# - nf-core/hisat2: This module provides the HISAT2 implementation.

# - nf-core/stringtie: This module provides the StringTie implementation.

# - nf-core/rsseqc: This module provides the RSeQC implementation.

# - nf-core/salmon_pseudo_aligner: This module provides the Salmon pseudo-aligner implementation.

# - nf-core/genomes: This module provides the genomes implementation.

# The pipeline is configured using the following parameters:

# - reads: The FASTQ file containing the reads.

# - reference: The FASTA file containing the reference genome.

# - outdir: The directory where the output files will be written.

# - salmon_index: The directory where the Salmon index will be written.

# - salmon_quant_libtype: The library preparation protocol (library type) used by Salmon quantification.

# - salmon_quant_args: Any additional arguments to be passed to Salmon quantification.

# - star_align_args: Any additional arguments to be passed to STAR alignment.

# - salmon_quant_args: Any additional arguments to be passed to Salmon quantification.

# - extra_star_align_args: Any additional arguments to be appended to the STAR align command.

# - extra_salmon_quant_args: Any additional arguments to be appended to the Salmon quant command.

# - seqBias: The seqBias parameter to be passed to Salmon quantification.

# - gcBias: The gcBias parameter to be passed to Salmon quantification.

# - cpus: The number of CPUs to use.

# - memory: The amount of memory to use.

# - cores: The number of cores to use.

# - threads: The number of threads to use.

# - cores_per_node: The number of cores per node.

# - nodes: The number of nodes to use.

# - max_job_size: The maximum job size to use.

# - min_job_size: The minimum job size to use.

# - max_runtime: The maximum runtime to use.

# - min_runtime: The minimum runtime to use.

# - max_walltime: The maximum walltime to use.

# - min_walltime: The minimum walltime to use.

# - loglevel: The log level to use.

# - debug: Whether to enable debug mode.

# - dryrun: Whether to run in dry run mode.

# - nohup: Whether to run in nohup mode.

# - output_format: The output format to use.

# - output_compression: The output compression to use.

# - output_directory: The output directory to use.

# - input_directory: The input directory to use.

# - working_directory: The working directory to use.

# - nfcore_config: The nf-core configuration file to use.

# - nfcore
