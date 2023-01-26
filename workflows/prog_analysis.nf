#!/bin/bash nextflow
params.outdir = 'results'

include { split_if_contains } from '../modules/functions.nf'
include {TCS}    from '../modules/evaluate_alignment.nf'

include {PROG_ALIGNER; PROG_ALIGNER_STRUCTURES; PROG_ALIGNER_EXPRESSO}       from '../modules/align.nf'
include {MMSEQS_PREP_DB; MMSEQS_SEARCH; TEMPLATE_FROM_DB_HITS; FETCH_STRUCTURES}        from '../modules/structures.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'
include {COLLECT_STRUCTURES_AF2} from '../subworkflows/collect_structures.nf'


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
    else if(align_method_string == "EXPRESSO"){

      // Collect the structures with
      MMSEQS_SEARCH(seqs_and_trees,target_db.collect())
      TEMPLATE_FROM_DB_HITS(MMSEQS_SEARCH.out.hits)
      FETCH_STRUCTURES(TEMPLATE_FROM_DB_HITS.out.filtered_hits)
      seqs_trees_structures = seqs_and_trees.combine(FETCH_STRUCTURES.out.fetched_structures, by: [0])

      PROG_ALIGNER_STRUCTURES (seqs_trees_structures, align_method)
      prog_alignment = PROG_ALIGNER_STRUCTURES.out.alignmentFile


    }else if(align_method_string == "3DCOFFEE"){


      COLLECT_STRUCTURES_AF2(seqs_and_trees.map{ it -> [it[0], it[1], "MAX", it[2]]}, structures)

      COLLECT_STRUCTURES_AF2.out.predicted_structures
      // Combine by family and tree
      seqs_and_trees
        .combine(COLLECT_STRUCTURES_AF2.out.predicted_structures, by: [0,1])
        .set{ seqs_and_trees_and_structures }

      PROG_ALIGNER_STRUCTURES (seqs_and_trees_and_structures, align_method)
      prog_alignment = PROG_ALIGNER_STRUCTURES.out.alignmentFile

      if (params.evaluate){
        EVALUATE_MSA( prog_alignment, refs_ch)
      }

    }else{

      PROG_ALIGNER (seqs_and_trees, align_method)
      prog_alignment = PROG_ALIGNER.out.alignmentFile
      TCS(prog_alignment)
      print("------------------------------------------------------")
      print(params.evaluate)
      prog_alignment.view()
      refs_ch.view()

      if (params.evaluate){
        print("here luisaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        prog_alignment.view()
        refs_ch.view()
        EVALUATE_MSA( prog_alignment, refs_ch)
      }

    }


    emit:
    alignment = prog_alignment
}
