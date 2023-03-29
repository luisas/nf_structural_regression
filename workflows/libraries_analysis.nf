#!/bin/bash nextflow
params.outdir = 'results'

include { split_if_contains } from '../modules/functions.nf'
include {TCS}    from '../modules/evaluate_alignment.nf'
include { TCS_WF } from '../modules/optimization.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'
include { MULTI_LIB_ALIGN } from '../subworkflows/library_align.nf'
include { ALIGN_WITH_LIBRARY } from '../modules/align.nf'

include { STRUCTURE_TO_3DI;   ENCODE_FASTA } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI; PROG_ALIGNER } from '../modules/align.nf'
include { MERGE_MAPPINGS } from '../modules/utils.nf'
include { FOLDSEEK_LIBRARY; SEQUENCE_LIBRARY; TMALIGN_LIBRARY;SAP_LIBRARY } from '../subworkflows/prep_libraries.nf'
include { MERGE_LIBRARIES } from '../modules/library.nf'
include { PRINTLIBRARY } from '../modules/library.nf'

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
    //seqs_and_trees_and_structures.view()
    fastas = seqs_and_trees.map { it -> [it[1], it[3]]}
    
    //seqs_and_trees_and_structures.view()
    if(library_method_string == "foldseek_only"){
        FOLDSEEK_LIBRARY(fastas, structures, matrix.collect())
        seqs_and_trees = seqs_and_trees.map{ it -> [ it[1], it[2], it[3], it[4]]}
        seqs_trees_libraries = seqs_and_trees.combine(FOLDSEEK_LIBRARY.out.library, by:0).view()
        ALIGN_WITH_LIBRARY(seqs_trees_libraries, library_method_string)
        if (params.evaluate){
             EVALUATE_MSA( ALIGN_WITH_LIBRARY.out.alignmentFile, refs)
        }
    }
    else if(library_method_string == "sequence_only" ){
        SEQUENCE_LIBRARY(fastas)
        seqs_and_trees = seqs_and_trees.map{ it -> [ it[1], it[2], it[3], it[4]]}
        seqs_trees_libraries = seqs_and_trees.combine(SEQUENCE_LIBRARY.out.library, by:0).view()
        ALIGN_WITH_LIBRARY(seqs_trees_libraries, library_method_string)
        if (params.evaluate){
             EVALUATE_MSA( ALIGN_WITH_LIBRARY.out.alignmentFile, refs)
        }
    }
    else if(library_method_string == "foldseek_sequence"){
        SEQUENCE_LIBRARY(fastas)
        FOLDSEEK_LIBRARY(fastas, structures, matrix.collect())
        MULTI_LIB_ALIGN(seqs_and_trees_and_structures, refs, library_method_string, FOLDSEEK_LIBRARY.out.library, SEQUENCE_LIBRARY.out.library)

    }else if(library_method_string == "tmalign_only"){
        TMALIGN_LIBRARY(seqs_and_trees_and_structures)
        seqs_and_trees = seqs_and_trees.map{ it -> [ it[1], it[2], it[3], it[4]]}
        seqs_trees_libraries = seqs_and_trees.combine(TMALIGN_LIBRARY.out.library, by:0).view()
        ALIGN_WITH_LIBRARY(seqs_trees_libraries, library_method_string)
        if (params.evaluate){
             EVALUATE_MSA( ALIGN_WITH_LIBRARY.out.alignmentFile, refs)
        }
    
    }else if(library_method_string == "sap_only"){
        SAP_LIBRARY(seqs_and_trees_and_structures)
        seqs_and_trees = seqs_and_trees.map{ it -> [ it[1], it[2], it[3], it[4]]}
        seqs_trees_libraries = seqs_and_trees.combine(SAP_LIBRARY.out.library, by:0).view()
        ALIGN_WITH_LIBRARY(seqs_trees_libraries, library_method_string)
        if (params.evaluate){
             EVALUATE_MSA( ALIGN_WITH_LIBRARY.out.alignmentFile, refs)
        }
    }
    else if(library_method_string == "tmalign_sequence"){
        SEQUENCE_LIBRARY(fastas)
        TMALIGN_LIBRARY(seqs_and_trees_and_structures)
        MULTI_LIB_ALIGN(seqs_and_trees, refs, library_method_string, TMALIGN_LIBRARY.out.library, SEQUENCE_LIBRARY.out.library)
    
    }else if(library_method_string == "sap_sequence"){
        SEQUENCE_LIBRARY(fastas)
        SAP_LIBRARY(seqs_and_trees_and_structures)
        MULTI_LIB_ALIGN(seqs_and_trees, refs, library_method_string, SAP_LIBRARY.out.library, SEQUENCE_LIBRARY.out.library)
    }

}
