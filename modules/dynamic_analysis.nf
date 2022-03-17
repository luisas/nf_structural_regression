#!/bin/bash nextflow

include {GENERATE_DYNAMIC_CONFIG; EXTRACT_SEQUENCES; ADD_PDB_HEADERS}      from './preprocess.nf'
include {REG_ALIGNER}       from './generateAlignment.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'
include {EVAL_ALIGNMENT}    from './modules_evaluateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'
include {HOMOPLASY}         from './modules_evaluateAlignment.nf'
include {METRICS}           from './modules_evaluateAlignment.nf'
include { RUN_COLABFOLD } from './localcolabfold.nf'


workflow DYNAMIC_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    tree_method
    bucket_size
    dynamicX

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

    // If AF2 predicted structures are to be used
    if (params.predict){
      // Extract sequences for prediction
      EXTRACT_SEQUENCES(seqs_and_trees, align_method, bucket_size, dynamicX, configFile, configValues, dynamicValues)
      // Predict with COLABFOLD
      RUN_COLABFOLD(EXTRACT_SEQUENCES.out.extractedSequences)
      // Prep PBD files for tcoffee
      ADD_PDB_HEADERS(RUN_COLABFOLD.out.af2_pdb)
      structures = ADD_PDB_HEADERS.out.pdb

    }else{
      structures = ""
    }

    // Then Align!
    DYNAMIC_ALIGNER (seqs_and_trees, align_method, bucket_size, dynamicX, configFile, configValues, dynamicValues, structures)



    /*
    *  Evaluate alignment
    */
    if (params.evaluate){
      refs_ch
        .cross (DYNAMIC_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT (alignment_and_ref)
      // Collect results in one CSV
      EVAL_ALIGNMENT.out.scores.map{ it -> "${it.baseName};${it.text}" }
                    .collectFile(name: "dynamic.scores.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")

    }

  emit:
  alignment = DYNAMIC_ALIGNER.out.alignmentFile
}
