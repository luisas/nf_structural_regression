
include {TCS;}    from '../modules/evaluate_alignment.nf'
include { split_if_contains } from '../modules/functions.nf'


workflow TCS_WF {

  take:
    alignments

  main:

    alignments.view()

    //EVAL_ALIGNMENT(alignment_and_ref)

}


// 1. Families with reference sequences and structures - channel (AF2 and experimental ? maybe af2 for now)
// 2. 3D coffee align
// 3. TCS

// if TCS < X:
// 4. remove sequence with lowest plddt
// 5. re-call 3D + TCS
