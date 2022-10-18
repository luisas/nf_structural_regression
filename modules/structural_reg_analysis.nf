#!/bin/bash nextflow

include {GENERATE_DYNAMIC_CONFIG}      from './preprocess.nf'
include {REG_ALIGNER}       from './generateAlignment.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'
include {EVAL_ALIGNMENT}    from './modules_evaluateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'
include {MMSEQS_PREP_DB; MMSEQS_SEARCH; TEMPLATE_FROM_DB_HITS; FETCH_STRUCTURES}        from './structures.nf'
include {STR_REG_ALIGNER}       from './generateAlignment.nf'


workflow STRUCTURAL_REG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
    target_db


  main:

    // 1. MMSEQS DATABASE PREPARE - nw only allowing to give the right one
    //MMSEQS_PREP_DB(target_db)

    // 2. Search against DB
    MMSEQS_SEARCH(seqs_and_trees,target_db.collect())

    // 3. Create template file
    TEMPLATE_FROM_DB_HITS(MMSEQS_SEARCH.out.hits)

    // 4. Download the PDBs
    FETCH_STRUCTURES(TEMPLATE_FROM_DB_HITS.out.filtered_hits)

    // Prep a channel with both sequence and structure infrmations
    seqs_trees_structures = seqs_and_trees.combine(FETCH_STRUCTURES.out.fetched_structures, by: [0])

    // 5. Align with regressive
    STR_REG_ALIGNER (seqs_and_trees, align_method, bucket_size)


}
