#!/bin/bash nextflow
params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process TREE_GENERATION {
    container 'edgano/tcoffee:pdb'
    tag "$tree_method on $id"
    storeDir "${params.outdir}/trees/${id}/${id}.${tree_method}"
    label "process_medium"

    input:
    tuple val(id), path(seqs)
    each tree_method

    output:
    tuple val (id), val(tree_method), path ("${id}.${tree_method}.dnd"), emit: trees
    path ".command.trace", emit: metricFile


    script:
    template "${path_templates}/tree/tree_${tree_method}.sh"
}
