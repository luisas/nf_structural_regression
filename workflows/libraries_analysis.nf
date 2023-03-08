#!/bin/bash nextflow
params.outdir = 'results'

include { split_if_contains } from '../modules/functions.nf'
include {TCS}    from '../modules/evaluate_alignment.nf'
include { TCS_WF } from '../modules/optimization.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'

include { ALIGN_WITH_LIBRARY } from '../modules/align.nf'

include { STRUCTURE_TO_3DI;   ENCODE_FASTA } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI } from '../modules/align.nf'
include {MERGE_MAPPINGS } from '../modules/utils.nf'
include {FOLDSEEK_LIBRARY; SEQUENCE_LIBRARY} from '../subworkflows/prep_libraries.nf'
include {MERGE_LIBRARIES} from '../modules/library.nf'
workflow LIBRARIES_ANALYSIS {
  take:
    seqs_and_trees
    refs
    library_method
    tree_method
    structures
    matrix

  main:

    library_method_string = library_method.toString().replace("[", "").replace("]", "")

    seqs_and_trees = seqs_and_trees.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1], it[2], it[3]]}
    seqs_and_trees_and_structures = seqs_and_trees.combine(structures, by: [0]).groupTuple(by:[1,2,3,4])
                                                    .map { it -> [ it[1], it[2], it[3], it[4], it[6]]}

    fastas = seqs_and_trees.map { it -> [it[0], it[3]]}

    if(library_method_string == "foldseek"){
        FOLDSEEK_LIBRARY(fastas, structures, matrix)
        // seqs_and_trees.combine(FOLDSEEK_LIBRARY.out.library, by:0).groupTuple(by:[1,2,3,4])
        //ALIGN_WITH_LIBRARY(seqs_and_trees,FOLDSEEK_LIBRARY.out.library, library_method_string)
        //EVAL 
    }
    else if(library_method_string == "foldseek_sequence"){
        SEQUENCE_LIBRARY(fastas)
        FOLDSEEK_LIBRARY(fastas, structures, matrix)
        
        // libraries = FOLDSEEK_LIBRARY.out.library.combine(SEQUENCE_LIBRARY.out.library, by:0)
        // MERGE_LIBRARIES(libraries,library_method_string)
        // ALIGN_WITH_LIBRARY(seqs_and_trees,MERGE_LIBRARIES.out.library, library_method)
        // EVAL
    }else if(library_method_string == "mtmalign_sequence"){
        //MTMALIGN_LIBRARY(seqs_and_trees, structures)
        //SEQUENCE_LIBRARY(seqs_and_trees, structures)
        //MERGE_LIBRARIES(MTMALIGN_LIBRARY.out.library, SEQUENCE_LIBRARY.out.library)
        //ALIGN_WITH_LIBRARY(seqs_and_trees,MERGE_LIBRARIES.out.library, library_method)
        // EVAL
    }





}
