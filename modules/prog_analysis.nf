#!/bin/bash nextflow
params.outdir = 'results'

include {EVAL_ALIGNMENT}      from './modules_evaluateAlignment.nf'
include {PROG_ALIGNER}       from './generateAlignment.nf'
workflow PROG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method

  main:
    PROG_ALIGNER (seqs_and_trees, align_method)

    if (params.evaluate){
      refs_ch
        .cross (PROG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }
      EVAL_ALIGNMENT (alignment_and_ref,align_method, bucket_size, dynamicX)
      EVAL_ALIGNMENT.out.scores.map{ it -> "${it.baseName};${it.text}" }
                    .collectFile(name: "progressive.scores.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")
    }

    emit:
    alignment = PROG_ALIGNER.out.alignmentFile
}
