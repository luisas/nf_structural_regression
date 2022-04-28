#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.dataset_dir="/users/cn/lsantus/"
//params.dataset_dir="/home/luisasantus/Desktop/crg_cluster"
params.af2_db_path = "${params.dataset_dir}/data/structural_regression/af2_structures"
params.fasta="$baseDir/a.fa"
Channel.fromPath("${params.af2_db_path}/**/*_alphafold.pdb").collectFile() { item ->
       [ "ids_done.txt", item.name + '\n' ]
   }.set{ids_done}

testfam = "test,seatoxin"
structures = Channel.fromPath("${params.af2_db_path}/colabfold_header/${testfam}/**/*.pdb")
structures_ch = structures.map { item -> [ item.getParent().getParent().baseName, item.baseName, item] }
structures_ch.view()
precomputed_structures_ids = structures_ch.collectFile() { item -> [ "ids_done.txt", item[1] + '\n' ]}.view()

//.map{it -> it.replace("id:", "")}

fastasplit = Channel.fromPath("${params.fasta}")
       .filter{ it.size()>0 }
       .splitFasta(record: [id: true])

fastasplit.view()
println(fastasplit.getClass())

fastasplit.map{it -> it.replace("id:", "")}.view()



workflow pipeline {
  //CHECK_CACHE("aa","test")
}


workflow {
  pipeline()
}
