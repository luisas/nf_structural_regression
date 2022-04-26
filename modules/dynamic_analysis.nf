#!/bin/bash nextflow

include {GENERATE_DYNAMIC_CONFIG; EXTRACT_SEQUENCES; ADD_PDB_HEADERS; CHECK_CACHE}      from './preprocess.nf'
include {REG_ALIGNER}       from './generateAlignment.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'
include {EVAL_ALIGNMENT}    from './modules_evaluateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'
include { RUN_COLABFOLD } from './localcolabfold.nf'


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

    /*
    *  Generate configuration for running dynamic t-coffee
    */
    GENERATE_DYNAMIC_CONFIG(dynamicMasterAln, dynamicSlaveAln, bucket_size, dynamicX)

    /*
    *  Align
    */

    // Test all configurations
    seqs_and_trees
        .combine(GENERATE_DYNAMIC_CONFIG.out.config)
        .set { seqs_and_trees_configs }



    if (params.predict){

      ids_done = structures.collectFile() { item -> [ "ids_done.txt", item[1] + '\n' ]}.collect()
      ids_done.view()
      // Extract sequences for which structures are needed
      EXTRACT_SEQUENCES(seqs_and_trees_configs)

      // 0. Check if structures have been produced
      CHECK_CACHE(EXTRACT_SEQUENCES.out.extractedSequences,ids_done)

      // 1A. PREDICT NEW STRUCTURES
      // Only runs if some sequences do not have a structure predicted
      RUN_COLABFOLD(CHECK_CACHE.out.seqToPredict.filter{ it[3].size()>0 })
      RUN_COLABFOLD.out.metricFile
                     .map{ it ->  "${it[0]}\n${it[1].text}\n${it[2].text}" }
                     .collectFile(name: "${workflow.runName}.trace", newLine: true, storeDir:"${params.af2_db_path}/colabfold/traces/")
      ADD_PDB_HEADERS(RUN_COLABFOLD.out.af2_pdb)


      // 1B. GET PRECOMPUTED STRUCTURES
      precomputed_structures = CHECK_CACHE.out.idsDone
                                          .map{ item -> [item.baseName.split("_")[0],
                                                         item.baseName.split("_")[1],
                                                         item.baseName.split("_")[2],
                                                         item.splitText().collect { it.trim()},
                                                         ]}
                                          .transpose()
                                          .map{ it -> [it[0], it[3], it[1], it[2]]}
                                          .join(structures, by: [0,1])
                                          .map{ it -> [it[0], it[2], it[3], it[4]]}

      // 2. Get a channel with both the newly and pre- computed structures
      all_structures = ADD_PDB_HEADERS.out.pdb.concat(precomputed_structures)
                                      .map{ it -> [it[0], it[1], it[3]] }
                                      .groupTuple(by: [0,1])

      //all_structures.view()
      seqs_and_trees
        .combine(all_structures, by: [0,1])
        .combine(EXTRACT_SEQUENCES.out.templates.map{ it -> [it[0], it[1], it[3], it[2], it[4], it[5], it[6], it[7]] }, by: [0,1])
        .set{ seqs_and_trees_and_structures }

    }else{
      empty = Channel.from('p','q').collect()
      structures = ""
    }


    //seqs_and_trees_and_structures.view()

    DYNAMIC_ALIGNER (seqs_and_trees_and_structures)


    /*
    *  Evaluate alignment
    */
    if (params.evaluate){
      refs_ch
        .cross (DYNAMIC_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT(alignment_and_ref)
      EVAL_ALIGNMENT.out.scores.collect().view()
      EVAL_ALIGNMENT.out.scores
                    .collectFile(name: "dynamic.scores.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")
    }



  emit:
  alignment = DYNAMIC_ALIGNER.out.alignmentFile
}
