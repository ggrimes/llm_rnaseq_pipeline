#!/usr/bin/env nextflow

// Load parameters from nextflow.config
params.reads = params.reads
params.index = params.index
params.out_dir = params.out_dir
params.num_cpus = params.num_cpus

// Define input channel for paired-end reads
Channel.fromFilePairs(params.reads)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { paired_reads }

// Define Salmon process for quantification
process salmon_quant {
    tag "$name"
    publishDir "${params.out_dir}/salmon", mode: 'copy'

    input:
    set val(name), file(reads) from paired_reads

    output:
    file "${name}_quant" into quant_ch

    script:
    """
    salmon quant -i ${params.index} -l A \
        -1 ${reads[0]} -2 ${reads[1]} \
        -p ${params.num_cpus} \
        -o ${name}_quant
    """
}

// Collect and merge quantification results
quant_ch.collect()
    .map { it -> [it.name, it / "quant.sf"] }
    .set { all_quant_results }

// Generate a combined gene expression matrix
process generate_matrix {
    publishDir "${params.out_dir}/expression_matrix", mode: 'copy'

    input:
    file(salmon_quant) from all_quant_results.flatten()

    output:
    file "expression_matrix.tsv"

    script:
    """
    tximport.R --input "${salmon_quant.join(',')}" \
        --output "expression_matrix.tsv"
    """
}
