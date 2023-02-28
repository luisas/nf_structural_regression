#!/bin/bash nextflow
params.outdir = 'results'

include { split_if_contains } from '../modules/functions.nf'
include {TCS}    from '../modules/evaluate_alignment.nf'
include { TCS_WF } from '../modules/optimization.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'

include { ALIGN_WITH_LIBRARY } from '../modules/align.nf'
include { GENERATE_LIBRARY } from '../modules/library.nf'

include { STRUCTURE_TO_3DI } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI } from '../modules/align.nf'

workflow FOLDSEEK_ANALYSIS {
  take:
    seqs_and_trees
    refs
    align_method
    tree_method
    structures

  main:

    align_method_string = align_method.toString().replace("[", "").replace("]", "")


    seqs_and_trees = seqs_and_trees.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1], it[2], it[3]]}
    seqs_and_trees_and_structures = seqs_and_trees.combine(structures, by: [0]).groupTuple(by:[1,2,3,4])
                                                    .map { it -> [ it[1], it[2], it[3], it[4], it[6]]}
    seqs_and_trees_and_structures.view()
    // tuple val(id), val(tree_method), path(seqs), path(guide_tree), path (structures)
    structures.view()

    if(align_method_string == "foldseek_library"){

      GENERATE_LIBRARY(seqs_and_trees_and_structures)

      ALIGN_WITH_LIBRARY(GENERATE_LIBRARY.out.library, align_method)
      alignment = ALIGN_WITH_LIBRARY.out.alignmentFile

    }else if(align_method_string == "foldseek_alphabet"){

      
      STRUCTURE_TO_3DI(structures.groupTuple(by:0))
      // ALIGN_WITH_3DI(FASTA_TO_3DI.out.fasta_3di)
      // DECODE_MSA(ALIGN_WITH_3DI.out.alignmentFile)
      // alignment = DECODE_MSA.out.alignmentFile
    }


    // if (params.evaluate){
    //   EVALUATE_MSA( alignment, refs)
    // }

  // emit: 
  //   alignmentFile = alignment

}
