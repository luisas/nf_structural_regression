#!/bin/bash nextflow
include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process GENERATE_DYNAMIC_CONFIG {
    container 'edgano/tcoffee:pdb'
    tag "Config 4 Dynamic"
    storeDir "${params.outdir}/dynamic_config/"
    label "small"

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
  tag "$align_method - $tree_method on $id; ${masterAln}-${masterSize}:${slaveAln}-${slaveSize}"
  //publishDir "${params.outdir}/seqs_extracted", pattern: '*'
  storeDir "${params.outdir}/seqs_extracted"

  input:
  tuple val(id), val(tree_method), path(seqs), path(guide_tree)
  val align_method
  each bucket_size
  each dynamicX
  path (dynamicConfig)
  tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
  val dynamicValues

  output:
  tuple val (id), path("*.fasta"), emit: extractedSequences

  script:
  template "${path_templates}/dynamic_align/dynamic_EXTRACT.sh"

}


process ADD_PDB_HEADERS{
  container 'edgano/tcoffee:pdb'
  tag "${fam_name}"
  storeDir "${params.outdir}/structures/colabfold_predictions_header/${fam_name}_colabfold"

  input:
  tuple val (fam_name), path (af2_pdb)

  output:
  tuple val(fam_name), path("*_header.pdb"), emit: pdb

  script:
  """
  for i in `find *.pdb`; do /tcoffee/t_coffee/src/extract_from_pdb -force -infile \$i > test.pdb; f="\$(basename -- \$i .pdb)"; mv test.pdb \${f}_header.pdb; done
  """
}
