
include {ADD_PDB_HEADERS; CHECK_CACHE}      from '../modules/preprocess.nf'
include { RUN_COLABFOLD } from '../modules/localcolabfold.nf'
include { split_if_contains } from '../modules/functions.nf'


workflow COLLECT_STRUCTURES_AF2 {

  take:
    seq_and_trees_mod
    structures

  main:

    ids_done = structures.ifEmpty(">--")
                               .collectFile() { item -> [ "ids_done.txt", item[1] + '\n' ]}
                               .collect()


    CHECK_CACHE(seq_and_trees_mod, ids_done)
    RUN_COLABFOLD(CHECK_CACHE.out.seqToPredict.filter{ it[3].size()>0 }.splitFasta( by: params.n_af2, file: true ),
                  params.model_type,
                  params.db )


    precomputed_structures = CHECK_CACHE.out.idsDone
                                        .map{ item -> [item.baseName.split("_")[0..-3].join("_"),
                                                       item.baseName.split("_")[-2],
                                                       item.baseName.split("_")[-1],
                                                       item.splitText().collect { it.trim()},
                                                       ]}
                                        .transpose()
                                        .map{ it -> [ split_if_contains(it[0], "-ref", 0), it[3].replace("/", "_"), it[1], it[2], it[0]]}
                                        .join(structures, by: [0,1])
                                        .map{ it -> [it[4], it[2], it[3], it[5]]}
    // 2. Get a channel with both the newly and pre- computed structures
    all_structures = RUN_COLABFOLD.out.af2_pdb.groupTuple(by:[0,1,2])
                                  .map{ it -> [it[0],it[1],it[2],it[3].flatten()]}
                                  .concat(precomputed_structures)
                                  .map{ it -> [it[0], it[1], it[3].toString().split("/")[-1].split("_alphafold")[0], it[3]] }


    // How it looks like: all structures -->  [fam,tree,sequence_id, sequence_id.pdb]
    ADD_PDB_HEADERS(all_structures.map{ it -> [ it[0], it[1], it[2], it[3], split_if_contains(it[0], "-ref", 0)]})



  emit:
    predicted_structures = ADD_PDB_HEADERS.out.pdb.groupTuple(by: [0,1])

}
