
nextflow.enable.dsl = 2

params.dataset_dir="/users/cn/lsantus/"
dataset = "homfam,extHomfam_v35-uniprot"

params.seqs ="${params.dataset_dir}/data/structural_regression/{${dataset}}/combinedSeqs/*.fa"
params.path_scripts = "$baseDir/bin"


params.outputdir = "${params.dataset_dir}/data/structural_regression/stats/"
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item.getParent().getParent().baseName, item] }


seqs_ch.view()


process CALC_SEQS_LENGTH{

  tag "${fam_name}"
  container 'luisas/python:bio3'
  storeDir "${params.outputdir}/seq_lengths/${dataset}"
  label "process_low"

  input:
  tuple val(fam_name), val(dataset), path(fasta)

  output:
  path ("${fam_name}_lengths.csv"),emit: seq_lengths

  script:
  template python3 "${params.path_scripts}/calc_seqlength.py"

}

process STATS_LENGTHS{

  container 'luisas/python:bio3'
  storeDir "${params.outputdir}/seq_lengths/"
  label "process_low"

  input:
  path(list_csvs)

  output:
  path ("summary_lengths.csv"),emit: stats

  script:
  template python3 "${params.path_scripts}/calc_seqlength_stats.py"

}



workflow pipeline {

    CALC_SEQS_LENGTH (seqs_ch)
    STATS_LENGTHS(CALC_SEQS_LENGTH.collect())

}

workflow {
  pipeline()
}
