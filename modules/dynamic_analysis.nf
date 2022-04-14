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
    dynamicX
    structures
    ids_done

  main:

    /*
    *  Generate configuration for running dynamic t-coffee
    */
    if(params.dynamicConfig){
      GENERATE_DYNAMIC_CONFIG(params.dynamicMasterAln, params.dynamicMasterSize, params.dynamicSlaveAln, params.dynamicSlaveSize)
      align_method="CONFIG"
      configFile = GENERATE_DYNAMIC_CONFIG.out.configFile
      configValues = GENERATE_DYNAMIC_CONFIG.out.configValues
      dynamicValues = "${params.dynamicMasterAln}.${params.dynamicMasterSize}_${params.dynamicSlaveAln}.${params.dynamicSlaveSize}"
    }else{
      align_method="DEFAULT"
      configFile = "/"
      configValues=["","","",""]
      dynamicValues = "DEFAULT"
    }


    /*
    *  Align
    */

    if (params.predict){
      // Extract sequences for which structures are needed
      EXTRACT_SEQUENCES(seqs_and_trees, align_method, bucket_size, dynamicX, configFile, configValues, dynamicValues)

      // 0. Check if structures have been produced
      ids_done.view()
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
        .combine(EXTRACT_SEQUENCES.out.templates.map{ it -> [it[0], it[1], it[3]] }, by: [0,1])
        .set{ seqs_and_trees_and_structures }

    }else{
      empty = Channel.from('p','q').collect()
      structures = ""
    }


    // Get everythin together ( sequences, references, trees, structures and parent sequences )
    // The latter is now just for making sure everything is fine, can be removed later!


    //seqs_and_trees_and_structures.view()
    DYNAMIC_ALIGNER (seqs_and_trees_and_structures, align_method, bucket_size, dynamicX, configFile, configValues, dynamicValues)



    /*
    *  Evaluate alignment
    */
    if (params.evaluate){
      refs_ch
        .cross (DYNAMIC_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      alignment_and_ref.view()
      EVAL_ALIGNMENT (alignment_and_ref)

      EVAL_ALIGNMENT.out.scores.view()
      EVAL_ALIGNMENT.out.scores.map{ it -> "${it.baseName};${it.text}" }
                    .collectFile(name: "dynamic.scores_${params.buckets}.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")
    }



  emit:
  alignment = DYNAMIC_ALIGNER.out.alignmentFile
}
