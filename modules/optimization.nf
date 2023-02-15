#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()




process TCS_WF2 {
  container 'luisas/tcoffee_python:2'
  tag "TCS on $id - ${test_alignment.baseName}"
  storeDir "${params.outdir}/evaluation/tcs_optimization/$id"
  label "process_low"

  input:
  tuple  val(id), val(align_method), file(fasta), file(structures), file (test_alignment)

  output:
  path ("*"), emit: tcs_wf_files
  path ("${id}_tcsoptimization_summary.csv"), emit: tcs_wf_summary

  script:
  """

  # Set starting variable names

  nfa=\$(grep ">" $fasta  | wc -l)
  NITER=1

  # Prep initial file names
  cp $fasta ${fasta.baseName}_\$NITER.fa
  cp ${test_alignment} ${test_alignment.baseName}_\$NITER.aln

  while [ \$nfa -gt 2 ]
  do

    FASTA=${fasta.baseName}_\$NITER.fa

    # 1. GET TOTAL TCS

    # Bad hack to circumvent t_coffee bug
    # Issue described already in: https://github.com/cbcrg/tcoffee/issues/3
    # Add an A in front of filename if the file begins with A
    filename=${test_alignment.baseName}_\$NITER.aln
    first_letter_filename=\${filename:0:1}
    if [ "\$first_letter_filename" == "A" ]; then input="A"\$filename; else input=\$filename;  fi

    t_coffee -infile \$input -evaluate -output=score_ascii -outfile "${test_alignment.baseName}_\$NITER.tcs"

    # 2. IF TCS_TOTAL > X: STOP, OTHERWISE: EXTRACT SEQ WITH MIN TCS

    TOTAL_TCS=\$(grep "SCORE" "${test_alignment.baseName}_\$NITER.tcs" | cut -d"=" -f2)
    # For now let's just do it systematically
    # Every time remove one and continue till there is something to align
    MIN_TCS=\$(python3 "${path_templates}/scripts/get_min_tcs.py" "${test_alignment.baseName}_\$NITER.tcs" "${test_alignment.baseName}_\$NITER")
    echo \$TOTAL_TCS","\$MIN_TCS","\$NITER >> ${id}_tcsoptimization_summary.csv

    # 3. Remove the sequence from the MSA
    NITER=\$(( \$NITER + 1 ))
    python3 "${path_templates}/scripts/remove_sequence.py" \$FASTA \$MIN_TCS > ${fasta.baseName}_\$NITER.fa

    # 4. Re-align

    # Prep templates
    for i in `awk 'sub(/^>/, "")' ${fasta.baseName}_\$NITER.fa`; do
        id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> template_list_\$NITER.txt
    done

    # Align
    t_coffee ${fasta.baseName}_\$NITER.fa -method TMalign_pair -template_file template_list_\$NITER.txt -output fasta_aln -outfile ${test_alignment.baseName}_\$NITER.aln

    # 5. Set variables for next iteration
    nfa=\$(grep ">"  ${fasta.baseName}_\$NITER.fa  | wc -l)


  done

  """
}




process TCS_WF {
  container 'luisas/tcoffee_python:2'
  tag "TCS on $id - ${test_alignment.baseName}"
  storeDir "${params.outdir}/evaluation/tcs_optimization_removestructure/$id/${test_alignment.baseName}"
  label "process_low"

  input:
  //tuple  val(id), val(align_method), file(fasta), file(structures), file (test_alignment)
  tuple  val(id), val(align_method), file(fasta), file(structures), file (test_alignment), file(ref_alignment)


  output:
  tuple path ("${id}*"), path("template*"), emit: tcs_wf_files
  path ("${test_alignment.baseName}_tcsoptimization_summary.csvll"), emit: tcs_wf_summary

  script:
  """

  # Set starting variable names
  NITER=1

  # Prep initial file names
  cp ${test_alignment} ${test_alignment.baseName}_\$NITER.aln

  for i in `awk 'sub(/^>/, "")' $fasta`; do
      id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> template_list_\$NITER.txt
  done

  nseqstemplate=\$(grep ">" template_list_\$NITER.txt  | wc -l)


  # -----------------   First entry is the input alignment -----------
  # Calc TCS and other metrics
  filename=${test_alignment.baseName}_\$NITER.aln
  first_letter_filename=\${filename:0:1}
  if [ "\$first_letter_filename" == "A" ]; then input="A"\$filename; else input=\$filename;  fi
  t_coffee -infile \$input -evaluate -output=score_ascii -outfile "${test_alignment.baseName}_\$NITER.tcs"
  TOTAL_TCS=\$(grep "SCORE" "${test_alignment.baseName}_\$NITER.tcs" | cut -d"=" -f2)

  TC=\$(t_coffee -other_pg aln_compare \
           -al1 ${ref_alignment} \
           -al2 ${test_alignment.baseName}_\$NITER.aln \
          -compare_mode tc \
          | grep -v "seq1" | grep -v '*' | \
          awk '{ print \$4}' ORS="\t" )



  SP=\$(t_coffee -other_pg aln_compare \
           -al1 ${ref_alignment} \
           -al2 ${test_alignment.baseName}_\$NITER.aln \
          -compare_mode sp \
          | grep -v "seq1" | grep -v '*' | \
          awk '{ print \$4}' ORS="\t" )

  echo \$TOTAL_TCS","-:-","\$NITER","\$SP","\$TC >> ${test_alignment.baseName}_tcsoptimization_summary.csv




  while [ \$nseqstemplate -gt 0 ]
  do

    TEMPLATE=template_list_\$NITER.txt

    # For now let's just do it systematically
    # Every time remove one and continue till there is something to align
    MIN_TCS=\$(python3 "${path_templates}/scripts/get_min_tcs.py" "${test_alignment.baseName}_\$NITER.tcs" "${test_alignment.baseName}_\$NITER" "\$TEMPLATE")
    echo "python3 ${path_templates}/scripts/get_min_tcs.py ${test_alignment.baseName}_\$NITER.tcs ${test_alignment.baseName}_\$NITER \$TEMPLATE" >> commands.sh
    MIN_SEQ=\$(echo \$MIN_TCS | cut -d':' -f1)
    #echo \$TOTAL_TCS","\$MIN_TCS","\$NITER >> ${id}_tcsoptimization_summary.csv

    NITER=\$(( \$NITER + 1 ))


    # ----------------------      S2         -----------------------------------
    # 4B. Re-align all the sequences, but do not use the structure of the worst sequence
    # Prep templates


    python3 ${path_templates}/scripts/remove_sequence_template.py \$TEMPLATE \$MIN_SEQ > template_list_\$NITER.txt

    t_coffee $fasta -method TMalign_pair -template_file template_list_\$NITER.txt -output fasta_aln -outfile ${test_alignment.baseName}_\$NITER.aln
    echo "t_coffee $fasta -method TMalign_pair -template_file template_list_\$NITER.txt -output fasta_aln -outfile ${test_alignment.baseName}_\$NITER.aln" >> commands.sh

    # --------------------------- Here also insert the evaluation

    # Bad hack to circumvent t_coffee bug
    # Issue described already in: https://github.com/cbcrg/tcoffee/issues/3
    # Add an A in front of filename if the file begins with A
    filename=${test_alignment.baseName}_\$NITER.aln
    first_letter_filename=\${filename:0:1}
    if [ "\$first_letter_filename" == "A" ]; then input="A"\$filename; else input=\$filename;  fi

    t_coffee -infile \$input -evaluate -output=score_ascii -outfile "${test_alignment.baseName}_\$NITER.tcs"

    # 2. IF TCS_TOTAL > X: STOP, OTHERWISE: EXTRACT SEQ WITH MIN TCS

    TOTAL_TCS=\$(grep "SCORE" "${test_alignment.baseName}_\$NITER.tcs" | cut -d"=" -f2)

    TC=\$(t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment.baseName}_\$NITER.aln \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" )



    SP=\$(t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment.baseName}_\$NITER.aln \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" )


    echo \$TOTAL_TCS","\$MIN_TCS","\$NITER","\$SP","\$TC >> ${test_alignment.baseName}_tcsoptimization_summary.csv


    # 5. Set variables for next iteration
    nseqstemplate=\$(grep ">"  template_list_\$NITER.txt | wc -l)


  done

  """
}
