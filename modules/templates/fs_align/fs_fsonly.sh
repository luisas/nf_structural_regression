t_coffee -seq $seqs \
         -method fs_pair \
         -template_file SELF_E_ \
         -template_dir_E_ $fs_dir \
         -outfile ${id}.fs_only.${tree_method}.aln \
         -output fasta_aln \
         -usetree $guide_tree \
         -fs_matrix $matrix 