#!/bin/bash nextflow

include {GENERATE_DYNAMIC_CONFIG}      from './preprocess.nf'
include {REG_ALIGNER}       from './generateAlignment.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'
include {EVAL_ALIGNMENT}    from './modules_evaluateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'

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
      refs_ch
        .cross (REG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }
      EVAL_ALIGNMENT (alignment_and_ref)
      EVAL_ALIGNMENT.out.scores.map{ it -> "${it.baseName};${it.text}" }
                    .collectFile(name: "regressive.scores_${params.buckets}.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")
    }

    emit:
    alignment = REG_ALIGNER.out.alignmentFile

}
