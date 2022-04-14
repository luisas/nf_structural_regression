#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.dataset_dir="/users/cn/lsantus/"
//params.dataset_dir="/home/luisasantus/Desktop/crg_cluster"
params.af2_db_path = "${params.dataset_dir}/data/structural_regression/af2_structures"
params.fasta="$baseDir/a.fa"
Channel.fromPath("${params.af2_db_path}/**/*_alphafold.pdb").collectFile() { item ->
       [ "ids_done.txt", item.name + '\n' ]
   }.set{ids_done}

testfam = "test"
structures = Channel.fromPath("${params.af2_db_path}/colabfold_header/${testfam}/**/*.pdb")
structures_ch = structures.map { item -> [ item.getParent().getParent().baseName, item.baseName, item] }
structures_ch.view()
precomputed_structures_ids = structures_ch.collectFile() { item -> [ "ids_done.txt", item[1] + '\n' ]}.view()


Channel.fromPath("${params.fasta}").filter{ it.size()>0 }.view()

params.a = "a"
process CHECK_CACHE{

  input:
  val(ids_done)
  val(fam)

  output:
  tuple val(fam), val(ids_done), emit: idsDone
  tuple val("${params.a}"), path("done.txt"), path(".command.trace"), emit: metricFile

  script:
  """
  echo ${ids_done} > done.txt
  echo ${params.a} > test.txt
  """

}



workflow pipeline {
  CHECK_CACHE("aa","test")
  CHECK_CACHE.out.metricFile
                 .map{ it ->  "${it[0]};${it[1].text};\n${it[2].text}" }
                 .collectFile(name: "trace", newLine: true, storeDir:"$baseDir/cachetesting/")


}


workflow {
  pipeline()
}
