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


process CHECK_CACHE{

  input:
  val(ids_done)
  val(fam)

  output:
  val(fam), val(ids_done), emit: idsDone

  script:
  """
  cat ${ids_done} > done.txt
  """

}



workflow pipeline {
  //CHECK_CACHE(ids_done,"test")
  params.list = "test.txt"
  Channel.fromPath("${params.list}")
        .map{ item -> [item.baseName, "A", item.splitText().collect { it.trim().replaceAll("_alphafold.pdb","") }]}
        .transpose()
        .set { ids_done }

  fam = "test"
  structures = Channel.fromPath("${params.af2_db_path}/colabfold_header/${fam}/**/*.pdb")
                      .map{ it -> [it.getParent().getParent().baseName, it.baseName, it]}.view()


  structures.join(ids_done, by: [0,1]).view()

  ids_done.map{ it -> [it[0], it[2]] }.groupTuple().view()
  // Channel.fromPath(params.list)
  //         .splitText()
  //         .map { file(params.af2_db_path + "/" +  it) }
  //         .set { file_list }
  //
  // file_list.view()

}


workflow {
  pipeline()
}
