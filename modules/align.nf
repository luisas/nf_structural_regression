#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process REG_ALIGNER {
    container 'edgano/tcoffee:pdb'
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
    //path "${id}.regressive.${bucket_size}.${align_method}.${tree_method}.bucket.log", emit: bucketLogFile


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
    tuple val (id), path ("${id}.progressive.${align_method}.${tree_method}.aln"), path ("*.lib"), emit: libraryFile
    path ".command.trace", emit: metricFile
    path "${id}.progressive.${align_method}.${tree_method}.tcs", emit: tcs_file
    //tuple val (id), val(align_method), file(seqs), path(structures), file("${id}.progressive.${align_method}.${tree_method}.aln"), emit: alignmentFiles

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
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), path(library)
    each align_method


    output:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), path(library), path("${id}.foldseek.${tree_method}.${align_method}.aln"), emit: alignment
    tuple val (id),path("${id}.foldseek.${tree_method}.${align_method}.aln"), emit: alignmentFile

    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -infile $seqs -lib $library -usetree=$guide_tree -outfile ${id}.foldseek.${tree_method}.${align_method}.aln -output fasta_aln
    """
}

process ALIGN_WITH_3DI {
    container 'luisas/foldseek_tcoffee:2'
    storeDir "${params.outdir}/alignments_libraries/$id/${id}.${library_method}.${tree_method}.${align_method}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(library), path(structures)
    each library_method

    output: 
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(library), path(structures), path("${id}.3di.${tree_method}.${align_method}.aln"), emit: alignmentFile

    script:
    """
    t_coffee -infile $seqs -lib $library -tree use=$tree_method -outfile ${id}.${library_method}.${tree_method}.${align_method}.aln -output fasta_aln
    """

}



process COMPACT_ALIGNER {
    container 'luisas/compact'
    tag "$align_method - $tree_method - $bucket_size on $id"
    storeDir "${params.outdir}/compact_benchmark/$id/${id}.regressive_comp_analysis.${bucket_size}.${align_method}.${tree_method}"
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
    //path "${id}.homoplasy", emit: homoplasyFile
    // path ".command.trace", emit: metricFile
    //path "${id}.regressive.${bucket_size}.${align_method}.${tree_method}.bucket.log", emit: bucketLogFile


    script:
    template "${path_templates}/compact_align/reg_${align_method}.sh"
}



process FS_ALIGNER {
    container 'luisas/fsmsa:3'
    tag "$id"
    storeDir "${params.outdir}/alignments/$id/${id}.regressive_fs_analysis.${library_method}.${tree_method}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), file(fs_dir)
    file(matrix)
    each(library_method)

    output:

    tuple val (id), path ("${id}.*.aln"), emit: alignmentFile
    tuple val (id), path ("${id}.*.aln"), path ("${id}.*.lib"), emit: libraryFile
    path ".command.trace", emit: metricFile


    script:
    template "${path_templates}/fs_align/fs_${library_method}.sh"
}

process FSREG_ALIGNER {
    container 'luisas/fsmsa:3'
    tag "$id"
    storeDir "${params.outdir}/regressive_foldseek/$id/${id}.regressive_foldseek_${params.targetDB}.${library_method}.${tree_method}.${bucket_size}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), file(templatedir)
    file(matrix)
    each(library_method)
    file(methodfile)
    each bucket_size


    output:
    tuple val (id), path ("${id}.*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile


    script:
    template "${path_templates}/regfs_align/regfs_template.sh"
}


process REG_ALIGNER_STR {
    container 'luisas/compact'
    tag "$align_method - $tree_method - $bucket_size on $id"
    storeDir "${params.outdir}/compact_benchmark/$id/${id}.regressive_comp_analysis.${bucket_size}.${align_method}.${tree_method}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree), path (structures)
    each align_method
    each bucket_size

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.*.aln"), emit: alignmentFile



    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"
}


process STRREG_ALIGNER {
    container 'luisas/fsmsa:3'
    tag "$id"
    storeDir "${params.outdir}/regressive_3d_test/$id/${id}.regressive_3d_${params.targetDB}.${library_method}.${tree_method}.${bucket_size}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), file(structures)
    each(align_method)
    file(methodfile)
    each bucket_size


    output:
    tuple val (id), path ("${id}.*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile


    script:
    template "${path_templates}/regfs_align/reg3d_template.sh"
}


process STRREG_ALIGNERTMALIGN {
    container 'luisas/fsmsa:3'
    tag "$id"
    storeDir "${params.outdir}/regressive_3dTmalign/$id/${id}.regressive_3dtmalign_${params.targetDB}.${align_method}.${tree_method}.${bucket_size}"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), file(structures)
    each(align_method)
    file(methodfile)
    each bucket_size


    output:
    tuple val (id), path ("${id}.*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile


    script:
    template "${path_templates}/progressive_align/prog_3DCOFFEE.sh"
}
