t_coffee -seq $seqs \
         -method fs_pair proba_pair \
         -template_file SELF_E_ \
         -template_dir_E_ $fs_dir \
         -outfile ${id}.fs_sequence.${tree_method}.aln \
         -output fasta_aln \
         -usetree $guide_tree \
         -fs_matrix $matrix \
         -out_lib ${id}.fs_sequence.${tree_method}.lib 