
#export NO_MAFFT_BINARIES=1
export MAFFT_BINARIES=''

#-thread ${task.cpus} \

t_coffee -reg -reg_method mafftginsi_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -n_core=1 \
         -reg_thread=1 \
         -outfile ${id}.regressive.${bucket_size}.${align_method}.${tree_method}.aln 2> tcoffee.stderr


