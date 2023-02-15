
include {EVAL_ALIGNMENT; EASEL_INFO; TCS; SIM}    from '../modules/evaluate_alignment.nf'
include { split_if_contains } from '../modules/functions.nf'


workflow EVALUATE_MSA {

  take:
    alignments
    references

  main:

    references
      .cross (alignments.map{ it -> [split_if_contains(it[0], "-ref", 0), it[1]]})
      .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
      .set { alignment_and_ref }

    alignment_and_ref.view()

    EVAL_ALIGNMENT(alignment_and_ref)

    if (params.evaluate_extended){
      EASEL_INFO(alignments)
      SIM(alignments)
    }



}
