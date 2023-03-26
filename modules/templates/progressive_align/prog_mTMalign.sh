# Prep selected list files
for i in `awk 'sub(/^>/, "")' ${seqs}`; do
    id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> 
done


do echo -e `grep "$i" "!{fam}"_selected_ref.template_list | awk '{print $3}'`'.pdb'; done > "!{fam}"_selected_ref.mtmalign
	mTM-align -i "!{fam}"_selected_ref.mtmalign -o "!{fam}"_selected_ref_mtmalign




 mTM-align mtmalign_structures_list.txt -o ${id}.progressive.${align_method}.${tree_method}.aln
