
{ time -p t_coffee -reg \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -thread ${task.cpus} \
         -expand \
         -outfile ${id}.regressive_comp_analysis.${bucket_size}.${align_method}.${tree_method}.aln 2> tcoffee.stderr ; } 2> time.txt

