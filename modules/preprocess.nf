#!/bin/bash nextflow
include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()

process GENERATE_DYNAMIC_CONFIG {
    container 'luisas/structural_regression'
    tag "Config 4 Dynamic"
    storeDir "${params.outdir}/dynamic_config/"
    label "small"

    input:
    val (masterAln)
    val (masterSize)
    val (slaveAln)
    val (slaveSize)

    output:
    path "${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config", emit: configFile
    tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize), emit: configValues

    script:
    """
    echo '${masterAln} ${masterSize}' > ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    echo '${slaveAln} ${slaveSize}' >> ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    """
}


process EXTRACT_SEQUENCES {
  container 'edgano/tcoffee:pdb'
  tag "$align_method - $tree_method on $id; ${masterAln}-${masterSize}:${slaveAln}-${slaveSize}"
  storeDir "${params.outdir}/seqs_extracted/${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}"

  input:
  tuple val(id), val(tree_method), path(seqs), path(guide_tree)
  val align_method
  each bucket_size
  each dynamicX
  path (dynamicConfig)
  tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
  val dynamicValues

  output:
  tuple val (id), val(tree_method),val(masterSize), path("*.fasta"), emit: extractedSequences
  path ".command.trace", emit: metricFile

  script:
  template "${path_templates}/dynamic_align/dynamic_EXTRACT.sh"

}


process ADD_PDB_HEADERS{
  container 'edgano/tcoffee:pdb'
  tag "${fam_name}"
  storeDir "${params.af2_db_path}/colabfold_header/${fam_name}/"

  input:
  tuple val (fam_name), val(tree_method),val(masterSize), path (af2_pdb)

  output:
  tuple val(fam_name), path("pdbs/*.pdb"), emit: pdb
  path("${fam_name}_plddt.eval"), emit: plddt

  script:
  """
  for i in `find *.pdb`; do /tcoffee/t_coffee/src/extract_from_pdb -force -infile \$i > test.pdb; f="\$(basename -- \$i .pdb)"; mv test.pdb \${f}_header.pdb; done
  # Store the
  for i in `find *_header.pdb`;do plddt=`awk '{print \$6"\t"\$11}' \$i | uniq | cut -f2 | awk '{ total += \$1 } END { print \$i total/(NR-1) }'`; echo \$i \$plddt >> plddt.eval; done
  sed 's/_alphafold_header.pdb//g' plddt.eval > ${fam_name}_plddt.eval

  mkdir pdbs
  for i in `find *_header.pdb`; do id_pdb=`echo \$i | sed 's._alphafold_header..g' | sed 's...g'`; mv \$i ./pdbs/\$id_pdb; done
  """
}


process CHECK_CACHE{

  container 'luisas/structural_regression'

  input:
  tuple val(fam_name),val(tree_method),val(dynamicMasterSize), path(fasta)
  val(ids_done)

  output:
  tuple val(fam_name),val(tree_method),val(dynamicMasterSize), path("toPredict.fasta"), emit: seqToPredict
  path("${fam_name}_${tree_method}_${dynamicMasterSize}.txt"), emit: idsDone

  shell:
  $/
  cat !{ids_done} > ids_done_full.txt

  # Extract all ids in the fasta file
  grep -o -E "^>.+" !{fasta} | tr -d ">" | sort | uniq > all_ids.txt

  # Extract all ids that already have a predicted structure
  sed 's/.pdb//g' ids_done_full.txt | sort | uniq > !{fam_name}_!{tree_method}_!{dynamicMasterSize}.txt

  # Check which ones need to be predicted
  diff --side-by-side --suppress-common-lines all_ids.txt !{fam_name}_!{tree_method}_!{dynamicMasterSize}.txt | cut -f1 > ids_to_predict.txt

  # Get fasta with only the sequences that were not found in af2_db
  awk 'NR==FNR{n[">"$$0];next} f{print f ORS $$0;f=""} $$0 in n{f=$$0}' ids_to_predict.txt !{fasta} > toPredict.fasta
  cat toPredict.fasta | tr -d '[:space:]' > toPredict.fasta
  /$


}
