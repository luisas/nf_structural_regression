#!/bin/bash nextflow
include { set_templates_path; split_if_contains  } from './functions.nf'
path_templates = set_templates_path()

process GENERATE_DYNAMIC_CONFIG {
    container 'luisas/structural_regression:7'
    tag "Config 4 Dynamic"
    storeDir "${params.outdir}/dynamic_config/"
    label "process_low"

    input:
    each masterAln
    each slaveAln

    output:
    tuple val(masterAln), val(slaveAln), path("${masterAln}_${slaveAln}.config"), emit: config


    script:
    """
    echo '${masterAln}' > ${masterAln}_${slaveAln}.config
    echo '${slaveAln}' >> ${masterAln}_${slaveAln}.config
    sed -i 's/#/\\n/g' ${masterAln}_${slaveAln}.config
    sed -i 's/:/\\t/g' ${masterAln}_${slaveAln}.config
    """
}


process EXTRACT_SEQUENCES {
  container 'luisas/structural_regression:7'
  tag "$id; ${id}.dynamic.${bucket_size}.${tree_method}"
  storeDir "${params.outdir}/seqs_extracted/${id}.dynamic.${bucket_size}.${tree_method}"
  label "process_low"


  input:
  tuple val(id), val(tree_method), path(seqs), path(guide_tree)
  each bucket_size

  output:
  tuple val (id), val(tree_method),val(bucket_size), path("${id}.${bucket_size}.${tree_method}.PARENTS.fasta"), emit: extractedSequences
  tuple val (id), val(tree_method),val(bucket_size), path("${id}.${bucket_size}.${tree_method}.templates.txt"), emit: templates
  path ".command.trace", emit: metricFile

  script:
  template "${path_templates}/dynamic_align/dynamic_EXTRACT.sh"

}


process ADD_PDB_HEADERS{
  container 'edgano/tcoffee:pdb'
  tag "${fam_name}"
  storeDir "${params.af2_db_path}/colabfold_header/${fam_name_short}/"
  label "process_low"

  input:
  tuple val (fam_name), val(tree_method), val(seq_id), path (af2_pdb), val(fam_name_short)

  output:
  tuple val(fam_name), val(tree_method), path("pdbs/${seq_id}.pdb"), emit: pdb
  path("plddts/${seq_id}_plddt.eval"), emit: plddt

  script:
  """
  # Add the headers
  mkdir pdbs
  /tcoffee/t_coffee/src/extract_from_pdb -force -infile $af2_pdb > pdbs/${seq_id}.pdb

  # Store the plddts summary
  mkdir plddts
  plddt=`awk -v var="${seq_id}" '{print \$6"\t"\$11}' pdbs/${seq_id}.pdb | uniq | cut -f2 | awk '{ total += \$1 } END { print var total/(NR-1) }'`; echo ${seq_id} \$plddt >> "plddts/${seq_id}_plddt.eval"
  """
}


process CHECK_CACHE{

  tag "${fam_name}"
  container 'luisas/structural_regression:7'
  label "process_low"


  input:
  tuple val(fam_name),val(tree_method),val(dynamicMasterSize), path(fasta)
  path(ids_done)

  output:
  tuple val(fam_name),val(tree_method),val(dynamicMasterSize), path("toPredict.fasta"), emit: seqToPredict
  path("${fam_name}_${tree_method}_${dynamicMasterSize}.txt"), emit: idsDone

  shell:
  $/
  cat !{ids_done} > ids_done_full.txt

  # Extract all ids in the fasta file
  grep -o -E "^>.+" !{fasta} | tr -d ">" | sort | uniq > all_ids.txt

  # Extract all ids that already have a predicted structure
  # sed 's/.pdb//g' ids_done_full.txt | sort | uniq > !{fam_name}_!{tree_method}_!{dynamicMasterSize}.txt

  # Check which ones need to be predicted and which ones already have a predicted structure
  touch !{fam_name}_!{tree_method}_!{dynamicMasterSize}.txt
  touch ids_to_predict.txt
  for id_fasta in `cat all_ids.txt`; do id_fasta_clean=`sed 's./._.g' <<< $$id_fasta`; if grep -Fqx $$id_fasta_clean ids_done_full.txt ; then echo $$id_fasta >> !{fam_name}_!{tree_method}_!{dynamicMasterSize}.txt; else echo $$id_fasta >> ids_to_predict.txt; fi;  done

  # Get fasta with only the sequences that were not found in af2_db
  awk 'NR==FNR{n[">"$$0];next} f{print f ORS $$0;f=""} $$0 in n{f=$$0}' ids_to_predict.txt !{fasta} > toPredict.fasta
  # cat toPredict.fasta | tr -d '[:space:]' > toPredict.fasta
  /$


}
