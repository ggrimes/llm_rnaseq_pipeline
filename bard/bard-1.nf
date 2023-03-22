
#!/usr/bin/env nextflow

# This pipeline uses Salmon to quantify RNA-Seq data.

# Define the reference genome and transcriptome.
params.reference = "/path/to/reference.fa"
params.transcriptome = "/path/to/transcriptome.fa"

# Define the output directory.
params.output = "/path/to/output"

# Run the pipeline.
nextflow run rnaseq-nf -with-docker
