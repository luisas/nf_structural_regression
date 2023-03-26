#!/bin/bash nextflow
include { FOLDSEEK_LIBRARY; SEQUENCE_LIBRARY; TMALIGN_LIBRARY;SAP_LIBRARY } from '../subworkflows/prep_libraries.nf'
include { split_if_contains } from '../modules/functions.nf'


workflow MSA_3DI_ANALYSIS {

  take:
    seqs_and_trees
    refs_ch
    library_method
    tree_method
    bucket_size
    structures
    matrix

  main:

    // Prep channels 
    print("Running the MSA_3DI_ANALYSIS workflow")
    seqs_and_trees = seqs_and_trees.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1], it[2], it[3]]}
    seqs_and_trees_and_structures = seqs_and_trees.combine(structures, by: [0]).groupTuple(by:[1,2,3,4])
                                                    .map { it -> [ it[1], it[2], it[3], it[4], it[6]]}
    fastas = seqs_and_trees.map { it -> [it[1], it[3]]}

    fastas.view()
    // Prep libraries 
    SEQUENCE_LIBRARY(fastas)
    FOLDSEEK_LIBRARY(fastas, structures, matrix.collect())
    //libraries = FOLDSEEK_LIBRARY.out.library.combine(SEQUENCE_LIBRARY.out.library, by:0)
    //MERGE_LIBRARIES(libraries,library_method_string)



}

