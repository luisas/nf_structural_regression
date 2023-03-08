
include {COMPACT_ALIGNER}  from '../modules/align.nf'
include {EVAL_ALIGNMENT; EASEL_INFO; TCS}    from '../modules/evaluate_alignment.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'


workflow COMPACT_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size


  main:

    COMPACT_ALIGNER (seqs_and_trees, align_method, bucket_size)

    if (params.evaluate){
      EVALUATE_MSA( COMPACT_ALIGNER.out.alignmentFile, refs_ch)
    }

    emit:
    alignment = COMPACT_ALIGNER.out.alignmentFile

}
