
p="\$(echo \${PWD} | sed 's#/#\\\\/#g')" # escape forward slashes
sed "s#FULLPATH#\${p}#g" fsproba.txt > method_tmp.txt
sed "s#TEMPDIR#$templatedir#g" method_tmp.txt > method.txt
t_coffee -reg -seq $seqs -nseq $bucket_size -method method.txt -quickcompact
