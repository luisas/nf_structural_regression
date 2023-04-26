
nextflow.enable.dsl = 2
dataset = "homfam"
params.dataset_dir="/users/cn/lsantus/"
params.dataset_dir="/users/cn/lsantus/"
params.seqs ="${params.dataset_dir}/data/structural_regression/${dataset}/combinedSeqs/*.fa"


refs = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item.getParent().getParent().baseName, item] }



process PREP_REF{

  tag "${fam_name}"
  container 'luisas/python:bio3'
  storeDir "${params.outputdir}/seq_lengths/${dataset}"
  label "process_low"

  input:
  tuple val(fam_name), val(dataset), path(fasta)

  output:
  path ("${dataset}_${fam_name}_lengths.csv"),emit: seq_lengths

  script:
  template "${params.path_scripts}/calc_seqlength.py"

}



workflow pipeline {
    CALC_SEQS_LENGTH (seqs_ch)
}

workflow {
  pipeline()
}
