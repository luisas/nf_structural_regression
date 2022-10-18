# Prep templates
for i in `awk 'sub(/^>/, "")' ${seqs}`; do
    id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> template_list.txt
done


t_coffee ${seqs} -method TMalign_pair -template_file "template_list.txt" -output fasta_aln -outfile ${id}.progressive.${align_method}.${tree_method}.aln
