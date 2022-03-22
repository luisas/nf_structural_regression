export blast_server_4_CLTCOFFEE=LOCAL
export NO_MAFFT_BINARIES=1
export VERBOSE_4_DYNAMIC=1
export DUMP_ALN_BUCKETS=1
export DUMP_SEQ_BUCKETS_ONLY=1


t_coffee -reg -reg_method dynamic_msa \
          -seq ${seqs} \
          -reg_tree ${guide_tree} \
          -reg_nseq ${bucket_size} \
          -dynamic ${dynamicX} \
          -dynamic_config ${dynamicConfig}

mv seqdump.1 ${id}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.PARENTS.fasta