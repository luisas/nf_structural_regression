#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process REG_ALIGNER {
    container 'luisas/structural_regression:16'
    tag "$align_method - $tree_method - $bucket_size on $id"
    storeDir "${params.outdir}/alignments/$id/${id}.regressive.${bucket_size}.${align_method}.${tree_method}"
    label 'process_medium'

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
    path "${id}.regressive.${bucket_size}.${align_method}.${tree_method}.bucket.log", emit: bucketLogFile


    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"
}


process DYNAMIC_ALIGNER {
    container 'luisas/structural_regression:17'
    tag "${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}"
    storeDir "${params.outdir}/alignments/$id/${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree), path (structures), path (extractedSequences),val(bucket_size), val(masterAln), val(slaveAln), path(dynamicConfig)
    each dynamicX

    output:
    tuple val (id), path("${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}.aln"), emit: alignmentFile
    path "*.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile
    path "${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}.bucket.log", emit: bucketLogFile

    script:
    template "${path_templates}/dynamic_align/dynamic_${masterAln}.sh"

}


process PROG_ALIGNER_STRUCTURES {
    container 'luisas/structural_regression:7'
    tag "$align_method - $tree_method on $id"
    storeDir "${params.outdir}/alignments/${id}/${id}.progressive.${align_method}_$db.${tree_method}"
    label 'process_small'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree), val(db), path(template), path (structures)
    each align_method


    output:
    tuple val (id), path ("${id}.progressive.${align_method}.${tree_method}.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}


process PROG_ALIGNER {
    container 'luisas/structural_regression:7'
    tag "$align_method - $tree_method on $id"
    storeDir "${params.outdir}/alignments/${id}/${id}.progressive.${align_method}.${tree_method}"
    label 'process_small'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method


    output:
    tuple val (id), path ("${id}.progressive.${align_method}.${tree_method}.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}



process STR_REG_ALIGNER {
    container 'luisas/structural_regression:17'
    tag "${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}"
    storeDir "${params.outdir}/alignments/$id/${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(db), path(structures)
    each align_method
    each bucket_size
    each dynamicX

    output:
    tuple val (id), path("${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}.aln"), emit: alignmentFile
    path "*.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile
    path "${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}.bucket.log", emit: bucketLogFile

    script:
    template "${path_templates}/dynamic_align/dynamic_${masterAln}.sh"

}
