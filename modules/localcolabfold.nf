#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RUN_COLABFOLD {
	tag "${fam_name}"
	storeDir "${params.af2_db_path}/colabfold/${fam_name}"

	input:
	tuple val(fam_name),val(tree_method),val(dynamicMasterSize), path(fasta)

	output:
	path ("*"), emit: all_af2
	tuple val(fam_name),val(tree_method),val(dynamicMasterSize),path ("*_alphafold.pdb"), emit: af2_pdb
	path ".command.trace", emit: metricFile

	script:
	"""
	colabfold_batch --amber --templates --num-recycle 3 ${fasta} \$PWD  ${params.cpu_flag}
	for i in `find *_relaxed_rank_1*.pdb`; do cp \$i `echo \$i | sed "s|_relaxed_rank_|\t|g" | cut -f1`"_alphafold.pdb"; done
	"""
}
