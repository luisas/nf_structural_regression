#! /bin/bash

declare compressFlag=" "

if $params.compressAZ ; then
    compressFlag=" -output fastaz_aln"
fi

{ time -p t_coffee -reg -reg_method t_coffee_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         \$compressFlag \
         -outfile ${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln 2> tcoffee.stderr ; } 2> time.txt

mv buckets.log ${id}.regressive.${bucket_size}.${align_method}.${tree_method}.bucket.log
