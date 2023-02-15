# Prep templates
for i in `awk 'sub(/^>/, "")' ${seqs}`; do
    id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> template_list.txt
done


t_coffee ${seqs} -method TMalign_pair -template_file "template_list.txt" -out_lib ${id}.progressive.${align_method}.${tree_method}.lib -output fasta_aln -outfile ${id}.progressive.${align_method}.${tree_method}.aln


# Calculate the TCS

filename=${id}
first_letter_filename=\${filename:0:1}
if [ "\$first_letter_filename" == "A" ]; then input="A"${id}.progressive.${align_method}.${tree_method}.aln; else input=${id}.progressive.${align_method}.${tree_method}.aln;  fi

t_coffee -infile \$input  -evaluate -lib ${id}.progressive.${align_method}.${tree_method}.lib -output score_ascii -outfile ${id}.progressive.${align_method}.${tree_method}.tcs
