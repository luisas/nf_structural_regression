#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process MMSEQS_PREP_DB {
    container 'soedinglab/mmseqs2'
    storeDir "${params.dbdir}/dbs/${params.target_db}"
    label 'process_small'
    tag "$db_id"

    input:
    tuple val(db_id), file(db)

    output:
    path db

    script:
    """
    mmseqs createindex $db/$db_id tmp
    """
}


process MMSEQS_SEARCH {
    container 'luisas/mmseqs2test'
    storeDir "${params.outdir}/structures/search_hits/mmseqs/$id/${id}.${db_id}"
    label 'process_small'
    tag "$id in $db_id"

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree)
    tuple val(db_id), file(db)

    output:
    tuple val(id), val(db_id), path("hits.m8"), emit: hits

    script:
    """
    mmseqs easy-search --min-seq-id 0.6 -c 0.7 --cov-mode 1 ${seqs} ${db}/${db_id} hits.m8 tmp
    """

}


process TEMPLATE_FROM_DB_HITS {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/search_hits/mmseqs/$id/${id}.${db_id}"
    label 'process_small'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits)

    output:
    tuple val(id), val(db_id), file("filtered_hits.m8"), file("template.txt"), file("ids_to_download.txt"), emit: filtered_hits

    script:
    """
    python3 "${path_templates}/scripts/filter_hits.py" ${hits} "filtered_hits.m8"
    """
}


// Currently bad, no caching but just or trying out things TODO --> fix
process FETCH_STRUCTURES {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/fetched/${db_id}/${id}/"
    label 'process_small'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download)

    output:
    tuple val(id), val(db_id), file(template), file("*.pdb"), emit: fetched_structures

    script:
    """
    for id in \$(cat $ids_to_download); do wget https://files.rcsb.org/download/\$id.pdb; done
    """
}
