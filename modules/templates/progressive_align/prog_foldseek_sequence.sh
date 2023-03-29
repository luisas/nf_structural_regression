t_coffee -seq $seqs \
         -method fs_pair proba_pair \
         -template_file SELF_E_ \
         -template_dir_E_ ${fs_mapping_dir} \
         -outfile ${id}_fs.aln \
         -output fasta_aln \
         -fs_matrix $matrix 