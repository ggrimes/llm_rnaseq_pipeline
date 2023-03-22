#!/usr/bin/env nextflow

# This pipeline uses Salmon to quantify RNA-Seq reads.

# Define the reference genome and transcriptome FASTA files.
params.reference = "/path/to/reference/genome.fa"
params.transcriptome = "/path/to/transcriptome.fa"

# Define the output directory.
params.outdir = "/path/to/output"

# Define the number of cores to use.
params.cpus = 8

# Define the output format for the results.
params.outputFormat = "tsv"

workflow {

  # Read the reads from the input FASTQ files.
  read_pairs_ch = channel.fromFilePairs( params.reads )

  # Index the reads with Salmon.
  index_ch = INDEX(read_pairs_ch)

  # Quantify the reads with Salmon.
  quant_ch = QUANT(index_ch)

  # Write the results to the output directory.
  write_results(quant_ch)
}

process INDEX {
  # Index the reads with Salmon.
  salmon index --threads $task.cpus -t $transcriptome -i ${read.index}
}

process QUANT {
  # Quantify the reads with Salmon.
  salmon quant --threads $task.cpus --libType=U -i ${read.index} - 1 ${read.first} - 2 ${read.second} -o ${read.id}
}

# Write the results to the output directory.
process write_results {
  # Write the results to a CSV file.
  writeResults(params.outdir, params.outputFormat)
}

# Write the results to a CSV file.
writeResults(params.outdir, params.outputFormat) {
  # Open a CSV file for writing.
  with CSVWriter(params.outdir + "results.tsv") {
    # Write the header row.
    write("Sample", "Reads", "Genes", "Isoforms", "Total")

    # Write the results for each sample.
    for (sample, reads, genes, isoforms, total) in results {
      write(sample, reads, genes, isoforms, total)
    }
  }
}
