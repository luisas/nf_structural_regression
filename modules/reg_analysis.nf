#!/bin/bash nextflow
//params.outdir = 'results_REG'
include {GENERATE_DYNAMIC_CONFIG}      from './preprocess.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'
include {EVAL_ALIGNMENT;EVAL_COMPRESS}    from './modules_evaluateAlignment.nf'


workflow DYNAMIC_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    tree_method
    bucket_size
    dynamicX

  main:
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

    DYNAMIC_ALIGNER (seqs_and_trees, align_method, bucket_size, dynamicX, configFile, configValues, dynamicValues)


    if (params.evaluate){
      refs_ch
        .cross (DYNAMIC_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      EVAL_ALIGNMENT ("dynamic", alignment_and_ref, DYNAMIC_ALIGNER.out.alignMethod, DYNAMIC_ALIGNER.out.treeMethod, DYNAMIC_ALIGNER.out.bucketSize)
      EVAL_ALIGNMENT.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.dynamic.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }

  emit:
  alignment = DYNAMIC_ALIGNER.out.alignmentFile
}
