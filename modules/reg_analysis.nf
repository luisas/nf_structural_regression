#!/bin/bash nextflow
//params.outdir = 'results_REG'
include {GENERATE_DYNAMIC_CONFIG}      from './preprocess.nf'
include {DYNAMIC_ALIGNER}             from './generateAlignment.nf'


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

  emit:
  alignment = DYNAMIC_ALIGNER.out.alignmentFile
}
