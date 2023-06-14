
nextflow.enable.dsl = 2
params.dataset = "benchfam"
params.dataset_dir_benchfam="/users/cn/lmansouri/clean_BENCHFAM"
params.outdir="/users/cn/lsantus/data/structural_regression/${params.dataset}"
params.testfam="PF16657,PF00004"

//params.seqs ="${params.dataset_dir_benchfam}/${params.testfam}/*_ref.aln"
params.seqs ="${params.dataset_dir_benchfam}/*/*_ref.aln"


ref_clustal = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.getParent().baseName, item] }

structures =  Channel.fromPath( "${params.dataset_dir_benchfam}/*/PDB/*.pdb" )
                     .map { item -> [ item.getParent().getParent().baseName,item.getParent().baseName, item] }
                     .groupTuple(by: [0,1])
templates =  Channel.fromPath( "${params.dataset_dir_benchfam}/*/*.template_list" )
                     .map { item -> [ item.getParent().baseName, item] }
                     .groupTuple(by: [0,1])



include {FOLDSEEK_CONVERT} from './workflows/foldseek.nf'
include { SAVE_MERGED_DIR } from './modules/utils.nf'

process PREP_REF{

  tag "${id}"
  container 'luisas/python:bio3'
  label "process_low"
  storeDir "${params.outdir}/refs"

  input:
  tuple val(id), path(aln)

  output:
  tuple path("${id}.ref"), path ("${id}-ref.fa"),emit: ref_fa

  script:
  """
  clustal_to_fa.py ${aln} ${id}.ref ${id}-ref.fa
  """

}


process PREP_STR{

  tag "${id}"
  container 'luisas/python:bio3'
  label "process_low"
  storeDir "${params.outdir}/structures/${id}"

  input:
  tuple val(id), val(db), path(structures), path(template)

  output:
  tuple val(id), val(db), path("${db}/*.pdb"),emit: structures_prep

  script:
  """
  mkdir -p ${db}

  rename_structures_with_template.py ${template} "make_links_tmp.sh" ${db}
  [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
  [ -f ./make_links.sh ] && bash ./make_links.sh
  [ -f ./make_links.sh ] && cat ./make_links.sh
  """

}



workflow pipeline {
  PREP_REF(ref_clustal)
  structures_and_templates = structures.combine(templates, by: [0])
  structures_renamed= PREP_STR(structures_and_templates)
  foldseek_db = FOLDSEEK_CONVERT(structures_renamed)
  foldseek_db = foldseek_db.map{  it -> [it[0].split("\\.")[0], it[1], it[2] ] }.groupTuple(by: [0,1])
  foldseek_db.view()
  SAVE_MERGED_DIR(foldseek_db, "foldseek")
}

workflow {
  pipeline()
}
