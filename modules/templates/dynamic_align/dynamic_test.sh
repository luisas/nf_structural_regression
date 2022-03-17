export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}
export NO_MAFFT_BINARIES=1
export VERBOSE_4_DYNAMIC=1



for i in `awk 'sub(/^>/, "")' ${seqs}`; do id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}"_alphafold_header.'pdb'; done  > template_list.txt

t_coffee ${seqs} \
         -method sap_pair,TMalign_pair \
         -template_file template_list.txt \
         -output fasta_aln \
         -outfile ${id}.dynamic_${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.with.${tree_method}.tree.aln \
         -thread 0
