
nextflow.enable.dsl = 2

seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item.getParent().getParent().baseName, item] }


process CALC_SEQS_LENGTH{

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


process SIM {
    container 'luisas/structural_regression:20'
    tag "SIM on $dataset - ${fasta.baseName}"
    storeDir "${params.outputdir}/sim_idscore/${dataset}"
    label "process_big"

    input:
    tuple  val(fam_name), val(dataset), path(fasta)

    output:
    tuple val(fam_name), val(dataset), path ("${fasta.baseName}.sim"), emit: sim

    script:
    """
    t_coffee -other_pg seq_reformat -in ${fasta} -output=sim_idscore > "${fasta.baseName}.sim"
    """
}


process STATS_LENGTHS{

  container 'luisas/python:bio3'
  publishDir "${params.outputdir}/seq_lengths/" , mode: 'copy'
  label "process_low"

  input:
  path(list_csvs)

  output:
  path ("summary_lengths.csv"),emit: stats

  script:
  template "${params.path_scripts}/calc_seqlength_stats.py"

}


process SIM_STATS {
    container 'luisas/structural_regression:20'
    tag "$fam_name"
    storeDir "${params.outputdir}/sim_idscore/${dataset}"
    label "process_small"

    input:
    tuple  val(fam_name), val(dataset), path(sim)

    output:
    path ("${sim.baseName}.sim_tot"), emit: sim_tot
    //path ("${sim.baseName}.sim_all"), emit: sim_all

    script:
    """
    echo "$fam_name" > tmp 
    echo "$dataset" >> tmp 
    grep ^TOT $sim | cut -f4 >> tmp
    tr '\n' ',' < tmp > ${sim.baseName}.sim_tot
    """
}


workflow pipeline {

    CALC_SEQS_LENGTH (seqs_ch)
    STATS_LENGTHS(CALC_SEQS_LENGTH.out.seq_lengths.toList())
    SIM(seqs_ch)
    SIM_STATS(SIM.out.sim)
    SIM_STATS.out.sim_tot.map{ it ->  "${it.text}" }
                    .collectFile(name: "${params.dataset}_similarities_summary.csv", newLine: true, storeDir:"${params.outdir}/sim_idscore/")  

}

workflow {
  pipeline()
}
