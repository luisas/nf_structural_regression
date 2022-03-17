#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process REG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'
    //publishDir "${params.outdir}/templates", pattern: '*.template_list'
    //publishDir "${params.outdir}/templates", pattern: '*.prf'

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
    path "time.txt", emit: timeFile
    //path "*.template_list", emit: templateFile
    //path "*.prf", emit: templateProfile

    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"
}

process PROG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.*.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile
    path "time.txt", emit: timeFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}

process DYNAMIC_ALIGNER {
    container 'luisas/structural_regression'
    //tag "$align_method - $tree_method on $id"
    tag "$align_method - $tree_method on $id; ${masterAln}-${masterSize}:${slaveAln}-${slaveSize}"
    storeDir "${params.outdir}/alignments/$id"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    val align_method
    each bucket_size
    each dynamicX
    path (dynamicConfig)
    tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
    val dynamicValues
    tuple val (fam_name), path (structures)

    output:
    tuple val (id), path("${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}.aln"), emit: alignmentFile
    //path "${id}.homoplasy", emit: homoplasyFile
    //path ".command.trace", emit: metricFile

    script:
    // Only use the MSA for the parent sequences since the one in the bottom is always the same
    template "${path_templates}/dynamic_align/dynamic_${align_method}_${masterAln}.sh"
    //template "${path_templates}/dynamic_align/dynamic_test.sh"
}
