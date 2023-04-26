# Prep templates
for i in `awk 'sub(/^>/, "")' ${seqs}`; do
    id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> template_list.txt
done

t_coffee -reg -reg_method 3dcoffee_msa \
        -seq ${seqs} \
        -reg_tree ${guide_tree} \
        -reg_nseq ${bucket_size} \
        -outfile mode1 -quickcompact \
        -method TMalign_pair  \
        -template_file "template_list.txt" \
        -outfile ${id}.regressive.${bucket_size}.${align_method}.${tree_method}.aln

