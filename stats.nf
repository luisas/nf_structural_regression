
nextflow.enable.dsl = 2

params.dataset_dir="/users/cn/lsantus/"
dataset = "homfam,extHomfam_v35-uniprot"
dataset = "homfam"
missing_fams = "ABC_tran,response_reg"
//missing_fams = "seatoxin,hip"

//params.seqs ="${params.dataset_dir}/data/structural_regression/{${dataset}}/combinedSeqs/*.fa"

params.seqs ="${params.dataset_dir}/data/structural_regression/${dataset}/combinedSeqs/*.fa"

params.path_scripts = "$baseDir/bin"
print("${params.dataset_dir}")

params.outputdir = "${params.dataset_dir}/data/structural_regression/stats/"
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item.getParent().getParent().baseName, item] }


//seqs_split = seqs_ch.splitFasta(by: 100_000, file: true).map { item -> [ item[0], item[1], item[2],item[2].baseName] }
//seqs_split.view()


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
    storeDir "${params.outdir}/sim/${dataset}"
    label "process_low"

    input:
    tuple  val(fam_name), val(dataset), path(fasta)

    output:
    tuple val(fam_name), val(dataset), path ("${fasta.baseName}.sim"), emit: sim

    script:
    """
    t_coffee -other_pg seq_reformat -in ${fasta} -output=sim > "${fasta.baseName}.sim"
    """
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
  template "${params.path_scripts}/calc_seqlength_stats.py"

}


process SIM_STATS {
    container 'luisas/structural_regression:20'
    tag "$fam_name"
    storeDir "${params.outdir}/sim/${dataset}"
    label "process_low"

    input:
    tuple  val(fam_name), val(dataset), path(sim)

    output:
    path ("${sim.baseName}.sim_tot"), emit: sim_tot

    script:
    """
    echo "$fam_name" > tmp 
    echo "$dataset" >> tmp 
    tail -n4 ${sim} | cut -f4 >> tmp
    tr '\n' ',' < tmp > ${sim.baseName}.sim_tot
    """
}





workflow pipeline {

    //CALC_SEQS_LENGTH (seqs_ch)
    //STATS_LENGTHS(CALC_SEQS_LENGTH.out.seq_lengths.toList())
    SIM(seqs_ch)
    SIM_STATS(SIM.out.sim)
    SIM_STATS.out.sim_tot.map{ it ->  "${it.text}" }
                    .collectFile(name: "${dataset}_similarities_summary.csv", newLine: false, storeDir:"${params.outdir}/sim/")  

}

workflow {
  pipeline()
}
