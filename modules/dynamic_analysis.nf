#!/bin/bash nextflow

include {GENERATE_DYNAMIC_CONFIG; EXTRACT_SEQUENCES; ADD_PDB_HEADERS; CHECK_CACHE}      from './preprocess.nf'
include {ADD_PDB_HEADERS as ADD_PDB_HEADERS_2}      from './preprocess.nf'

include {REG_ALIGNER}       from './generateAlignment.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'
include {EVAL_ALIGNMENT}    from './modules_evaluateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'
include { RUN_COLABFOLD } from './localcolabfold.nf'
include { split_if_contains } from './functions.nf'


workflow DYNAMIC_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    tree_method
    bucket_size
    dynamicMasterAln
    dynamicSlaveAln
    dynamicX
    structures

  main:

    if (params.predict){

      println("here")
      // Extract sequences for which structures are needed
      EXTRACT_SEQUENCES(seqs_and_trees, bucket_size)

      // Add pdb header to available structures
      // The third value is empty because it corresponds to the masterSize. Redesign needed in the future.ifEmpty(">--")

      ids_done = structures.ifEmpty(">--")
                           .collectFile() { item -> [ "ids_done.txt", item[1] + '\n' ]}
                           .collect()

      // 0. Check if structures have been produced
      CHECK_CACHE(EXTRACT_SEQUENCES.out.extractedSequences, ids_done)

      // 1A. PREDICT NEW STRUCTURES
      // Only runs if some sequences do not have a structure predicted
      // Split Fasta and run colabfold on chunks
      RUN_COLABFOLD(CHECK_CACHE.out.seqToPredict.filter{ it[3].size()>0 }.splitFasta( by: params.n_af2, file: true ))



      // 1B. GET PRECOMPUTED STRUCTURES
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

      // Combine by family and tree
      seqs_and_trees
        .combine(ADD_PDB_HEADERS.out.pdb.groupTuple(by: [0,1]), by: [0,1])
        .combine(EXTRACT_SEQUENCES.out.templates.map{ it -> [it[0], it[1], it[3], it[2]] }, by: [0,1])
        .set{ seqs_and_trees_and_structures }

    }else{
       // TO BE TESTED!!
      structures = Channel.empty().ifEmpty("--")
                           .collectFile() { item -> ["DUMMY_FILE"]}
                           .collect()
      seqs_and_trees
        .cross(all_structures)
        .combine(EXTRACT_SEQUENCES.out.templates.map{ it -> [it[0], it[1], it[3]] }, by: [0,1])
        .set{ seqs_and_trees_and_structures }

    }


    /*
    *  Generate configuration for running dynamic t-coffee
    */
    GENERATE_DYNAMIC_CONFIG(dynamicMasterAln, dynamicSlaveAln)

    seqs_and_trees_and_structures
    .combine(GENERATE_DYNAMIC_CONFIG.out.config)
    .set { seqs_and_trees_and_structures_configs }


    DYNAMIC_ALIGNER (seqs_and_trees_and_structures_configs, dynamicX)


    /*
    *  Evaluate alignment
    */
    if (params.evaluate){


      refs_ch
        .cross (DYNAMIC_ALIGNER.out.alignmentFile.map{ it -> [split_if_contains(it[0], "-ref", 0), it[1]]})
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }


      EVAL_ALIGNMENT(alignment_and_ref)
      EASEL_INFO (alignment_and_ref)

      EVAL_ALIGNMENT.out.scores.collect()
      EVAL_ALIGNMENT.out.scores
                    .collectFile(name: "dynamic.scores.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")
    }



  emit:
  alignment = DYNAMIC_ALIGNER.out.alignmentFile
}
