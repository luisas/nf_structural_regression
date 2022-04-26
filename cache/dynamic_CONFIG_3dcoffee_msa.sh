export VERBOSE_4_DYNAMIC=1
export DUMP_ALN_BUCKETS=1




# -------- Create template file
# -------- Only extract parent sequences for the template!
# Extract parent sequences
#parent_sequences=`grep ">" ${extractedSequences} | sed 's.>..g'`

for i in `awk 'sub(/^>/, "")' ${seqs}`; do
  if grep -Fx ">\$i" ${extractedSequences}; then
    #id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}"_alphafold_header.'pdb' >> template_list.txt
    id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\$PWD/\${id_pdb}" >> template_list.txt
  else
    echo "-"
  fi
done


# -------- Run alignment
t_coffee -reg -reg_method dynamic_msa \
          -seq ${seqs} \
          -reg_tree ${guide_tree} \
          -reg_nseq ${bucket_size} \
          -dynamic ${dynamicX} \
          -dynamic_config ${dynamicConfig} \
          -output fasta_aln \
          -thread 0 \
          -reg_homoplasy \
          -template_file template_list.txt \
          -method TMalign \
          -outfile ${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}.aln



mv *.homoplasy ${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${masterSize}.${slaveAln}.${slaveSize}.${tree_method}.homoplasy
#-reg_homoplasy \
