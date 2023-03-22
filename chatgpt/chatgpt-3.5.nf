// Define the input directory containing the FASTQ files
params.input_dir = "./fastq"

// Define the output directory
params.output_dir = "./salmon"

// Define the reference transcriptome
params.transcriptome = "reference_transcriptome.fa"

// Define the number of threads to use for Salmon
params.salmon_threads = 8

// Define the Salmon index output directory
params.salmon_index_dir = "./salmon_index"

// Define the Nextflow process for indexing the transcriptome with Salmon
process salmon_index {
  input:
  file(transcriptome) from params.transcriptome

  output:
  dir(params.salmon_index_dir)

  script:
  """
  salmon index -t ${transcriptome} -i ${params.salmon_index_dir}
  """
}

// Define the Nextflow process for quantifying gene expression with Salmon
process salmon_quant {
  input:
  file(fastq) from channel.fromPath("${params.input_dir}/*.fastq.gz")

  output:
  file("${params.output_dir}/${fastq.baseName}.quant.sf")

  script:
  """
  salmon quant -i ${params.salmon_index_dir} -l A -p ${params.salmon_threads} \
    -r ${fastq} -o ${params.output_dir}/${fastq.baseName}_quant
  """
}

// Define the Nextflow workflow
workflow {
  // Run the Salmon index process first
  salmon_index()

  // Then run the Salmon quant process on each FASTQ file
  fastq_files = Channel.fromPath("${params.input_dir}/*.fastq.gz")
  salmon_quant(fastq_files)
}
