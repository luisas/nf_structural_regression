export blast_server_4_CLTCOFFEE=LOCAL
export NO_MAFFT_BINARIES=1
export VERBOSE_4_DYNAMIC=1
export DUMP_ALN_BUCKETS=1
export DUMP_SEQ_BUCKETS_ONLY=1


t_coffee -reg -reg_method dynamic_msa \
          -seq ${seqs} \
          -reg_tree ${guide_tree} \
          -reg_nseq ${bucket_size} \
          -thread ${task.cpus}

mv seqdump.1 ${id}.${bucket_size}.${tree_method}.PARENTS.fasta


for i in `awk 'sub(/^>/, "")' ${id}.${bucket_size}.${tree_method}.PARENTS.fasta`; do
    id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> ${id}.${bucket_size}.${tree_method}.templates.txt
done
