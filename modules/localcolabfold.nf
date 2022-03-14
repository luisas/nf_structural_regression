#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RUN_COLABFOLD {
	tag "${fam_name}"
	//publishDir "${params.outdir}/${fam_name}_colabfold", mode: 'copy'
	storeDir "${params.outdir}/${fam_name}_colabfold"

	input:
	tuple val(fam_name), path(fasta)

	output:
	path ("*"), emit: all_af2
	path ("*_alphafold.pdb"), emit: af2_pdb

	script:
	if (params.cpu_flag == true) {cpu_flag = '--cpu'}
	else {cpu_flag = ' '}
	"""
	colabfold_batch --amber --templates --num-recycle 3 ${fasta} \$PWD ${cpu_flag}
	for i in `find *_relaxed_rank_1*.pdb`; do cp \$i `echo \$i | sed "s|_relaxed_rank_|\t|g" | cut -f1`"_alphafold.pdb"; done
	"""
}
