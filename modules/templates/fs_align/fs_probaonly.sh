t_coffee -seq $seqs \
         -method proba_pair \
         -template_file SELF_E_ \
         -template_dir_E_ $fs_dir \
         -outfile ${id}.proba_only.${tree_method}.aln \
         -output fasta_aln \
         -usetree $guide_tree \
         -fs_matrix $matrix \
         -out_lib ${id}.proba_only.${tree_method}.lib