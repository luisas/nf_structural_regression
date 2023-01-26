#!/bin/bash nextflow

include {REG_ALIGNER}  from '../modules/align.nf'
include {EVAL_ALIGNMENT; EASEL_INFO; TCS}    from '../modules/evaluate_alignment.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'


workflow REG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size


  main:

    REG_ALIGNER (seqs_and_trees, align_method, bucket_size)

    if (params.evaluate){

      EVALUATE_MSA( REG_ALIGNER.out.alignmentFile, refs_ch)

    }

    emit:
    alignment = REG_ALIGNER.out.alignmentFile

}
