#! /bin/bash

declare compressFlag=" "

if $params.compressAZ ; then
    compressFlag=" -output fastaz_aln"
fi

t_coffee -reg -reg_method probcons_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_homoplasy \
         -thread ${task.cpus} \
         \$compressFlag \
         -outfile ${id}.regressive.${bucket_size}.${align_method}.${tree_method}.aln 2> tcoffee.stderr
