#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process REG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $bucket_size on $id"
    storeDir "${params.outdir}/alignments/$id/${id}.regressive.${bucket_size}.${align_method}.${tree_method}"

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree)
    each align_method
    each bucket_size

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.*.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"
}


process DYNAMIC_ALIGNER {
    container 'luisas/structural_regression:7'
    tag "$align_method - $tree_method on $id; ${masterAln}-${masterSize}:${slaveAln}-${slaveSize}"
    storeDir "${params.outdir}/alignments/$id/${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}"
    label 'process_small'
    afterScript 'sleep 10'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree), path (structures), path (extractedSequences)
    val align_method
    each bucket_size
    each dynamicX
    path (dynamicConfig)
    tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
    val dynamicValues

    output:
    tuple val (id), path("${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}.aln"), emit: alignmentFile
    path "*.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/dynamic_align/dynamic_${align_method}_${masterAln}.sh"

}


process PROG_ALIGNER {
    container 'luisas/structural_regression:7'
    tag "$align_method - $tree_method on $id"
    storeDir "${params.outdir}/alignments/${id}/${id}.progressive.${align_method}.${tree_method}"

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method

    output:
    tuple val (id), path ("${id}.progressive.*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}
