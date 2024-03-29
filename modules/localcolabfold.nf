process RUN_COLABFOLD {
	tag "${fam_name}"
	publishDir "${params.af2_db_path}/colabfold/${fam_name}", mode: 'copy', overwrite: true
	container "/users/cn/abaltzis/sing_cache/athbaltzis-colabfold_proteinfold-v0.9.img"
	input:
	tuple val(fam_name),val(tree_method),val(dynamicMasterSize), path(fasta)
	val   model_type
	path db

	output:
	path ("*"), emit: all_af2
	tuple val(fam_name),val(tree_method),val(dynamicMasterSize),path ("*_alphafold.pdb"), emit: af2_pdb
	tuple val(fam_name),path("*_hash_colabfold_trace.txt"), path(".command.trace"), emit: metricFile

	script:
	"""
	grep ">" ${fasta} > colabfold_trace.txt
	echo "#cpu flag" ${params.cpu_flag} >> colabfold_trace.txt
	echo "RUNNING TIME" >> colabfold_trace.txt
	hash=`sha1sum  colabfold_trace.txt |  head -c 40`

	{ time -p colabfold_batch --amber --templates --num-recycle 3 \
				 --data ${db}/params/${params.model_type} \
				 --model-type ${model_type} \
				 ${fasta} \
				 \$PWD  ${params.cpu_flag} 2> std.err ; } 2> time.txt
	cat time.txt >> colabfold_trace.txt
	mv colabfold_trace.txt "\${hash}"_hash_colabfold_trace.txt

	for i in `find *_relaxed_rank_1*.pdb`; do cp \$i `echo \$i | sed "s|_relaxed_rank_|\t|g" | cut -f1`"_alphafold.pdb"; done
	"""
}
