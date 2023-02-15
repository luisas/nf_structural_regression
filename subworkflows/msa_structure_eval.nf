
include {EXTRACT_MSA_PAIRS;GET_DISTANCES; GET_SCORE; }    from '../modules/evaluate_alignment.nf'
include { split_if_contains } from '../modules/functions.nf'


workflow MSA_IRMSD {

  take:
    references
    structures

  main:


    //refs_ch.combine( structures_ch, by: 0 ).groupTuple(by : [0,1]).view()
    // 1. Extract pairs
    EXTRACT_MSA_PAIRS(references)
    pairs_with_structures = EXTRACT_MSA_PAIRS.out.pairs.transpose().combine(structures, by:0)
    pairs_with_structures.view()
    // 4. Compute iRMSD
    // takes pairwise alignment --> output matrix with iRMSD per residue
    //GET_DISTANCES(pairs_with_structures, params.evaluation_metric, local_radius)
    //GET_SCORE(GET_DISTANCES.out.distances)
    // I want the final score in a format of ID1 ID2 POS SCORE
    //GET_SCORE.out.score.groupTuple(by:[0,1,2]).view()
    //MERGE_SCORES( GET_SCORE.out.score.groupTuple(by:[0,1,2]))


}
