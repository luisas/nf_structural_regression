p="\$(echo \${PWD} | sed 's#/#\\\\/#g')" # escape forward slashes
sed "s#FULLPATH#\${p}#g" $methodfile > method.txt
sed "s#/#_#g" $seqs > input.fa
sed "s#/#_#g" $guide_tree > tree.dnd

for i in `awk 'sub(/^>/, "")' input.fa`; do
    idpdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\$idpdb" >> template_list.txt
done


t_coffee -reg -seq input.fa -nseq $bucket_size -reg_tree tree.dnd -method method.txt -quickcompact -outfile ${id}.3d_regressive.${bucket_size}.${align_method}.${tree_method}.aln


