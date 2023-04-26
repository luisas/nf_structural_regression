#!/bin/bash nextflow
include { FOLDSEEK_LIBRARY; SEQUENCE_LIBRARY; TMALIGN_LIBRARY;SAP_LIBRARY } from '../subworkflows/prep_libraries.nf'
include { split_if_contains } from '../modules/functions.nf'
include { COMPACT_ALIGNER; FS_ALIGNER; FSREG_ALIGNER; STRREG_ALIGNER} from '../modules/align.nf'
include { STRUCTURE_TO_3DI;   ENCODE_FASTA } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI } from '../modules/align.nf'
include { MERGE_MAPPINGS } from '../modules/utils.nf'
include { PREP_FS_SEQS } from '../modules/encoding.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'

workflow STRUCTURAL_ANALYSIS {

  take:
    seqs_and_trees
    refs
    library_method
    tree_method
    bucket_size
    structures
    foldseek_structures
    matrix
    methodfile
    methodfile_3d

  main:

    // Prep channels 
    seqs_and_trees = seqs_and_trees.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1], it[2], it[3]]}
    seqs_and_trees_and_structures = seqs_and_trees.combine(structures, by: [0]).groupTuple(by:[1,2,3,4])
                                                    .map { it -> [ it[1], it[2], it[3], it[4], it[6]]}

    seqs_and_trees_and_foldseek = seqs_and_trees.combine(foldseek_structures, by: [0]).groupTuple(by:[1,2,3,4])
                                                    .map { it -> [ it[1], it[2], it[3], it[4], it[6]]}

    //seqs_and_trees_and_structures.view()                                                

    // 1. ALIGN WITH FOLDSEEK
    FSREG_ALIGNER(seqs_and_trees_and_foldseek, matrix.collect(), library_method, methodfile.collect(), bucket_size )

    // 2. ALIGN WITH STRUCTURES FULL 
    //STRREG_ALIGNER (seqs_and_trees_and_structures, "3dcoffee",methodfile_3d.collect(), bucket_size)

    // // EVALUATE 
    if (params.evaluate){
      EVALUATE_MSA( FSREG_ALIGNER.out.alignmentFile, refs)
    }

}

