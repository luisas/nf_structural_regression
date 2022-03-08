#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process DYNAMIC_ALIGNER {
    container 'edgano/tcoffee:pdb'
    //tag "$align_method - $tree_method on $id"
    tag "$align_method - $tree_method on $id; ${masterAln}-${masterSize}:${slaveAln}-${slaveSize}"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    val align_method
    each bucket_size
    each dynamicX
    path (dynamicConfig)
    tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
    val dynamicValues

    output:
    val "${dynamicValues}_${params.db}", emit: alignMethod
    val tree_method, emit: treeMethod
    val "${bucket_size}_${dynamicX}", emit: bucketSize
    tuple val (id), path("*.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    // template can be default or using a config file
    template "${path_templates}/dynamic_align/dynamic_${align_method}.sh"
}
