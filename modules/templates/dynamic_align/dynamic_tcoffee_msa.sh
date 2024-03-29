export VERBOSE_4_DYNAMIC=1
export method_4_CLTCOFFEE=TMalign_pair
export template_file_4_CLTCOFFEE="template_list.txt"
export QUIET_ENV=1
export LOG_4_DYNAMIC=1

# -------- Create template file
# -------- Only extract parent sequences for the template!
# Extract parent sequences
#parent_sequences=`grep ">" ${extractedSequences} | sed 's.>..g'`
mv *.templates.txt template_list.txt

# -------- Run alignment
t_coffee -reg -reg_method dynamic_msa \
          -seq ${seqs} \
          -reg_tree ${guide_tree} \
          -reg_nseq ${bucket_size} \
          -dynamic ${dynamicX} \
          -dynamic_config ${dynamicConfig} \
          -output fasta_aln \
          -thread ${task.cpus} \
          -reg_homoplasy \
          -outfile ${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}.aln



mv *.homoplasy ${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}.homoplasy
mv buckets.log ${id}.dynamic.${bucket_size}.dynamicX.${dynamicX}.${masterAln}.${slaveAln}.${tree_method}.bucket.log
