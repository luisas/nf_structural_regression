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
    container 'luisas/structural_regression:20'
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
    container 'luisas/structural_regression:20'
    tag "$align_method - $tree_method on $id"
    storeDir "${params.outdir}/alignments/${id}/${id}.progressive.${align_method}.${tree_method}"
    label 'process_small'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree), path (structures)
    each align_method


    output:
    tuple val (id), path ("${id}.progressive.${align_method}.${tree_method}.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile
    path "${id}.progressive.${align_method}.${tree_method}.tcs", emit: tcs_file
    tuple val (id), val(align_method), file(seqs), file(structures), file("${id}.progressive.${align_method}.${tree_method}.aln"), emit: alignmentFiles

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}


process PROG_ALIGNER {
    container 'luisas/structural_regression:20'
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

process PROG_ALIGNER_EXPRESSO {
    container 'luisas/expresso:latest'
    tag "$align_method - $tree_method on $id"
    storeDir "${params.outdir}/alignments/${id}/${id}.progressive.${align_method}.${tree_method}"
    label 'process_small'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method
    tuple val(db_id), path(db), path(index)


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

process ALIGN_WITH_LIBRARY{
    container 'luisas/structural_regression:20'
    storeDir "${params.outdir}/alignments_foldseek/$id/${id}.foldseek.${tree_method}.${align_method}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(library), path(structures)
    each align_method


    output:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(library), path(structures), path("${id}.foldseek.${tree_method}.${align_method}.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -infile $seqs -lib $library -use tree=$tree_method -outfile ${id}.foldseek.${tree_method}.${align_method}.aln -output fasta_aln
    """
}

process ALIGN_WITH_3DI {
    container 'luisas/structural_regression:20'
    storeDir "${params.outdir}/alignments_foldseek/$id/${id}.3di.${tree_method}.${align_method}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(library), path(structures)
    each align_method

    output: 
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(library), path(structures), path("${id}.3di.${tree_method}.${align_method}.aln"), emit: alignmentFile

    script:
    """
    t_coffee -infile $seqs -lib $library -outfile ${id}.3di.${tree_method}.${align_method}.aln -output fasta_aln
    """

}