

nextflow.enable.dsl = 2

include { RUN_COLABFOLD } from './modules/localcolabfold'


params.mode == "colabfold"

params.outdir = "$baseDir/results_storedir"
smallfam="seatoxin,hip"
//params.dataset_dir="/users/cn/lsantus/"
//params.dataset_dir="/home/luisasantus/Desktop/crg_cluster"
//params.seqs ="${params.dataset_dir}/data/structural_regression/homfam/combinedSeqs/{${smallfam}}.fa"
params.seqs = "${baseDir}/data/*.fasta"
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

seqs_ch.view()

workflow pipeline {
    if(params.mode == "colabfold") {
          RUN_COLABFOLD(seqs_ch)
    }
}



workflow {
  pipeline()
}

println "Done!"
