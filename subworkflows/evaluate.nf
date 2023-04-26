
include {EVAL_ALIGNMENT; EASEL_INFO; TCS; SIM; GAPS_PROGRESSIVE; EVAL_IRMSD}    from '../modules/evaluate_alignment.nf'
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

    EVAL_ALIGNMENT(alignment_and_ref)

    if (params.evaluate_extended){
      GAPS_PROGRESSIVE(alignments)
    }
}


workflow EVALUATE_MSA_STRUCTURAL {

  take:
    alignments
    structures

  main:

    alignment_and_structures = alignments.map{ it -> [split_if_contains(it[0], "-ref", 0), it[1]]}.combine(structures, by: [0]).groupTuple(by:[0,1])
                                                    .map { it -> [ it[0], it[1], it[2]]}
    
    EVAL_IRMSD(alignment_and_structures)


}