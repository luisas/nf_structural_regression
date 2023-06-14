#!/bin/bash nextflow
include { FOLDSEEK_LIBRARY; SEQUENCE_LIBRARY; TMALIGN_LIBRARY;SAP_LIBRARY } from '../subworkflows/prep_libraries.nf'
include { split_if_contains } from '../modules/functions.nf'
include { COMPACT_ALIGNER; FS_ALIGNER; FSREG_ALIGNER; STRREG_ALIGNER; STRREG_ALIGNERTMALIGN; PROG_ALIGNER_STRUCTURES} from '../modules/align.nf'
include { STRUCTURE_TO_3DI; } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI } from '../modules/align.nf'
include { MERGE_MAPPINGS } from '../modules/utils.nf'
include { PREP_FS_SEQS } from '../modules/encoding.nf'
include { EVALUATE_MSA ; EVALUATE_MSA_STRUCTURAL; EVALUATE_MSA_LIBRARIES;  } from '../subworkflows/evaluate.nf'
include { EVALUATE_MSA as EVALUATE_MSA_2 ; EVALUATE_MSA_STRUCTURAL as EVALUATE_MSA_STRUCTURAL_2 ; EVALUATE_MSA_LIBRARIES as EVALUATE_MSA_LIBRARIES_2;  } from '../subworkflows/evaluate.nf'

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
                                     
    // 1. ALIGN WITH FOLDSEEK
    //FSREG_ALIGNER(seqs_and_trees_and_foldseek, matrix.collect(), library_method, methodfile.collect(), bucket_size )
    FS_ALIGNER(seqs_and_trees_and_foldseek, matrix.collect(), library_method)

    // 2. ALIGN WITH STRUCTURES FULL 
    PROG_ALIGNER_STRUCTURES (seqs_and_trees_and_structures, "3DCOFFEE")
    //STRREG_ALIGNERTMALIGN (seqs_and_trees_and_structures, "3dcoffee",methodfile_3d.collect(), bucket_size)

    FS_ALIGNER.out.alignmentFile.view()

    // EVALUATE 
    if (params.evaluate){
      EVALUATE_MSA( FS_ALIGNER.out.alignmentFile, refs)
      EVALUATE_MSA_STRUCTURAL ( FS_ALIGNER.out.alignmentFile, structures)
      EVALUATE_MSA_LIBRARIES ( FS_ALIGNER.out.libraryFile )

      EVALUATE_MSA_2( PROG_ALIGNER_STRUCTURES.out.alignmentFile, refs)
      EVALUATE_MSA_STRUCTURAL_2 ( PROG_ALIGNER_STRUCTURES.out.alignmentFile, structures)
      EVALUATE_MSA_LIBRARIES_2 ( PROG_ALIGNER_STRUCTURES.out.libraryFile )
    }

    //final_alignment_file = FS_ALIGNER.out.alignmentFile.filter { it -> it[1] =~ "/*fs_sequence*/" }
    //final_alignment_file.view()
    //all_libraries = FS_ALIGNER.out.libraryFile.map{ it  -> [ it[0], it[2]]}.groupTuple(by: [0]).flatten().collect().view()
    //final_alignment_file.combine(all_libraries, by:0).view()
}

