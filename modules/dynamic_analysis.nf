#!/bin/bash nextflow

include {GENERATE_DYNAMIC_CONFIG}      from './preprocess.nf'
include {REG_ALIGNER}       from './generateAlignment.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'
include {EVAL_ALIGNMENT;EVAL_COMPRESS}    from './modules_evaluateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'
include {HOMOPLASY}         from './modules_evaluateAlignment.nf'
include {METRICS}           from './modules_evaluateAlignment.nf'
workflow REG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size

  main:
    //seqs_and_trees.view()
    REG_ALIGNER (seqs_and_trees, align_method, bucket_size)

    if (params.evaluate){
      refs_ch
        .cross (REG_ALIGNER.out.alignmentFile)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
        .set { alignment_and_ref }

      if (params.compressAZ){
      EVAL_COMPRESS ("regressive", alignment_and_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
      EVAL_COMPRESS.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_COMPRESS.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_COMPRESS.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      }else{
      EVAL_ALIGNMENT ("regressive", alignment_and_ref, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize)
      EVAL_ALIGNMENT.out.tcScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.spScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      EVAL_ALIGNMENT.out.colScore
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
      }

    }
    if (params.homoplasy){
      HOMOPLASY("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.homoplasyFile)
      HOMOPLASY.out.homoFiles
                    .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text};${it[10].text}" }
                    .collectFile(name: "${workflow.runName}.regressive.homo.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }

    def metrics_regressive = params.metrics? METRICS("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize, REG_ALIGNER.out.metricFile, REG_ALIGNER.out.timeFile) : Channel.empty()
    if (params.metrics) {
        metrics_regressive.metricFiles
                          .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text};${it[6].text};${it[7].text};${it[8].text};${it[9].text}" }
                          .collectFile(name: "${workflow.runName}.regressive.metrics.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }

    def easel_info = params.easel? EASEL_INFO ("regressive", REG_ALIGNER.out.alignmentFile, REG_ALIGNER.out.alignMethod, REG_ALIGNER.out.treeMethod, REG_ALIGNER.out.bucketSize) : Channel.empty()
    if (params.easel) {
        easel_info.easelFiles
                  .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[6].text};${it[7].text}" }
                  .collectFile(name: "${workflow.runName}.regressive.easel.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    }

    emit:
    alignment = REG_ALIGNER.out.alignmentFile
    metrics = metrics_regressive
    easel = easel_info

}

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
    DYNAMIC_ALIGNER (seqs_and_trees, align_method, bucket_size, dynamicX, configFile, configValues, dynamicValues)


    /*
    *  Evaluate alignment
    */
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
