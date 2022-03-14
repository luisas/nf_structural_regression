#!/bin/bash nextflow
params.outdir = 'results'
include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process GENERATE_DYNAMIC_CONFIG {
    container 'edgano/tcoffee:pdb'
    tag "Config 4 Dynamic"

    input:
    val (masterAln)
    val (masterSize)
    val (slaveAln)
    val (slaveSize)

    output:
      path "${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config", emit: configFile
      tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize), emit: configValues

    script:
    """
    echo '${masterAln} ${masterSize}' > ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    echo '${slaveAln} ${slaveSize}' >> ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    """
}


process EXTRACT_SEQUENCES {
  container 'edgano/tcoffee:pdb'
  //tag "$align_method - $tree_method on $id"
  tag "$align_method - $tree_method on $id; ${masterAln}-${masterSize}:${slaveAln}-${slaveSize}"
  publishDir "${params.outdir}/seqs_extracted", pattern: '*'

  input:
  tuple val(id), val(tree_method), path(seqs), path(guide_tree)
  val align_method
  each bucket_size
  each dynamicX
  path (dynamicConfig)
  tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
  val dynamicValues

  output:
  val "${dynamicValues}_${params.db}", emit: alignMethod
  val tree_method, emit: treeMethod
  val "${bucket_size}_${dynamicX}", emit: bucketSize
  tuple val (id), path("*"), emit: extractedSequences

  script:
  template "${path_templates}/dynamic_align/dynamic_EXTRACT.sh"

}
