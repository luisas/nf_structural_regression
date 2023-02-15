#!/bin/bash nextflow
params.outdir = 'results'

include { split_if_contains } from '../modules/functions.nf'
include {TCS}    from '../modules/evaluate_alignment.nf'

include {PROG_ALIGNER; PROG_ALIGNER_STRUCTURES; PROG_ALIGNER_EXPRESSO}       from '../modules/align.nf'
include {MMSEQS_PREP_DB; MMSEQS_SEARCH; TEMPLATE_FROM_DB_HITS; FETCH_STRUCTURES}        from '../modules/structures.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'
include {COLLECT_STRUCTURES_AF2} from '../subworkflows/collect_structures.nf'
include { TCS_WF } from '../modules/optimization.nf'


workflow PROG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    structures
    target_db
    channeled_dbs

  main:

    align_method_string = align_method.toString().replace("[", "").replace("]", "")

    if(align_method_string == "EXPRESSOREAL"){

      PROG_ALIGNER_EXPRESSO (seqs_and_trees, align_method, channeled_dbs)
      prog_alignment = PROG_ALIGNER_EXPRESSO.out.alignmentFile

    }
    else if(align_method_string == "3DCOFFEEEXPERIMENTAL"){

      seqs_and_trees = seqs_and_trees.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1], it[2], it[3]]}
      seqs_and_trees_and_structures = seqs_and_trees.combine(structures, by: [0]).groupTuple(by:[1,2,3,4])
                                                    .map { it -> [ it[1], it[2], it[3], it[4], it[6]]}
      seqs_and_trees_and_structures.view()
      // tuple val(id), val(tree_method), path(seqs), path(guide_tree), path (structures)


      PROG_ALIGNER_STRUCTURES (seqs_and_trees_and_structures, align_method)
      prog_alignment = PROG_ALIGNER_STRUCTURES.out.alignmentFile
      prog_alignment_files = PROG_ALIGNER_STRUCTURES.out.alignmentFiles

      refs_ch
        .cross (prog_alignment_files.map{ it -> [split_if_contains(it[0], "-ref", 0), it[1], it[2], it[3], it[4]]})
        .map { it -> [ it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[0][1] ] }
        .set{prog_alignment_and_ref}

      TCS(prog_alignment)
      TCS_WF(prog_alignment_and_ref)
      if (params.evaluate){
        EVALUATE_MSA( prog_alignment, refs_ch)
      }


    } else if(align_method_string == "3DCOFFEE"){


      COLLECT_STRUCTURES_AF2(seqs_and_trees.map{ it -> [it[0], it[1], "MAX", it[2]]}, structures)

      seqs_and_trees
        .combine(COLLECT_STRUCTURES_AF2.out.predicted_structures, by: [0,1])
        .set{ seqs_and_trees_and_structures }

      PROG_ALIGNER_STRUCTURES (seqs_and_trees_and_structures, align_method)
      prog_alignment = PROG_ALIGNER_STRUCTURES.out.alignmentFile
      prog_alignment_files = PROG_ALIGNER_STRUCTURES.out.alignmentFiles

      //prog_alignment_files.view()
      //refs_ch.view()

      refs_ch
        .cross (prog_alignment_files.map{ it -> [split_if_contains(it[0], "-ref", 0), it[1], it[2], it[3], it[4]]})
        .map { it -> [ it[1][0], it[1][1], it[1][2], it[1][3], it[1][4], it[0][1] ] }
        .set{prog_alignment_and_ref}

      TCS(prog_alignment)
      TCS_WF(prog_alignment_and_ref)

      if (params.evaluate){
        EVALUATE_MSA( prog_alignment, refs_ch)
      }

    }else{

      PROG_ALIGNER (seqs_and_trees, align_method)
      prog_alignment = PROG_ALIGNER.out.alignmentFile
      TCS(prog_alignment)

      if (params.evaluate){
        EVALUATE_MSA( prog_alignment, refs_ch)
      }

    }


    emit:
    alignment = prog_alignment
}
