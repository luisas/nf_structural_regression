#!/bin/bash nextflow
params.outdir = 'results'

include {EVAL_ALIGNMENT}      from './modules_evaluateAlignment.nf'
include {PROG_ALIGNER}       from './generateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'



workflow PROG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method

  main:
    PROG_ALIGNER (seqs_and_trees, align_method)

    if (params.evaluate){

      alignments = PROG_ALIGNER.out.alignmentFile.map{ it -> [ it[0].replaceAll("-ref", ""), it[1] ] }

      refs_ch
        .cross (alignments)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }.view()
        .set { alignment_and_ref }
      EVAL_ALIGNMENT (alignment_and_ref)
      EASEL_INFO (alignment_and_ref)

      EVAL_ALIGNMENT.out.scores
                    .collectFile(name: "progressive.scores.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")
    }

    emit:
    alignment = PROG_ALIGNER.out.alignmentFile
}
