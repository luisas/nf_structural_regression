
p="\$(echo \${PWD} | sed 's#/#\\\\/#g')" # escape forward slashes
sed "s#FULLPATH#\${p}#g" $methodfile > method_tmp.txt
sed "s#TEMPDIR#$templatedir#g" method_tmp.txt > method.txt
sed "s#/#_#g" $seqs > input.fa
sed "s#/#_#g" $guide_tree > tree.dnd

t_coffee -reg -seq input.fa -nseq $bucket_size -reg_tree tree.dnd -method method.txt -quickcompact -outfile ${id}.foldeek_regressive.${bucket_size}.${library_method}.${tree_method}.aln


