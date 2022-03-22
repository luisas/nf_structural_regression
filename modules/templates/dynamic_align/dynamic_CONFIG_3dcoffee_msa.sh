export blast_server_4_CLTCOFFEE=LOCAL
export protein_db_4_CLTCOFFEE=${params.database_path}
export NO_MAFFT_BINARIES=1
export VERBOSE_4_DYNAMIC=1
export DUMP_ALN_BUCKETS=1




# -------- Create template file
# -------- Only extract parent sequences for the template!
for i in `awk 'sub(/^>/, "")' ${seqs}`; do id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}"_alphafold_header.'pdb'; done  > template_list.txt

# -------- Run alignment
t_coffee -reg -reg_method dynamic_msa \
          -seq ${seqs} \
          -reg_tree ${guide_tree} \
          -reg_nseq ${bucket_size} \
          -dynamic ${dynamicX} \
          -reg_homoplasy \
          -dynamic_config ${dynamicConfig} \
          -template_file template_list.txt \
          -output fasta_aln \
          -outfile ${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}.aln

mv *.homoplasy ${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}.homoplasy
