#!/bin/bash nextflow
params.outdir = 'results'

process EVAL_ALIGNMENT {
    container 'luisas/structural_regression:7'
    tag "EVAL_ALIGNMENT on $id"
    storeDir "${params.outdir}/evaluation/score/${test_alignment.baseName}"
    label "process_low"

    input:
    tuple  val(id), file (test_alignment), file (ref_alignment)

    output:
    path ("${test_alignment.baseName}.scores"), emit: scores


    script:
    """
    ## Sum-of-Pairs Score ##
    t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            >> "scores.txt"

    ## Total Column Score ##
    t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            >> "scores.txt"

    ## Column Score ##
        t_coffee -other_pg aln_compare \
              -al1 ${ref_alignment} \
              -al2 ${test_alignment} \
             -compare_mode column \
             | grep -v "seq1" | grep -v '*' | \
             awk '{ print \$4}' ORS="\t" \
             >> "scores.txt"

    cat scores.txt | tr -s '[:blank:]' ';'  >  "${test_alignment.baseName}.scores"
    sed -i 's/^/${test_alignment.baseName};/' "${test_alignment.baseName}.scores"
    """
}


process TCS {
    container 'luisas/structural_regression:20'
    tag "TCS on $id - ${test_alignment.baseName}"
    storeDir "${params.outdir}/evaluation/tcs/"
    label "process_low"

    input:
    tuple  val(id), file (test_alignment)

    output:
    path ("${test_alignment.baseName}.tcs"), emit: tcs_score

    script:
    """

    # Bad hack to circumvent t_coffee bug
    # Issue described already in: https://github.com/cbcrg/tcoffee/issues/3
    # Add an A in front of filename if the file begins with A
    filename=${test_alignment.fileName}
    first_letter_filename=\${filename:0:1}
    if [ "\$first_letter_filename" == "A" ]; then input="A"\$filename; else input=\$filename;  fi

    t_coffee -infile \$input -evaluate -output=score_ascii -outfile "${test_alignment.baseName}.tcs"

    """
}

process TCS_ON_LIB {
    container 'luisas/structural_regression:20'
    tag "TCS on $id - ${test_alignment.baseName}"
    storeDir "${params.outdir}/evaluation/tcs/"
    label "process_low"

    input:
    tuple  val(id), file (test_alignment), file(lib)

    output:
    path ("${test_alignment.baseName}.tcs"), emit: tcs_score

    script:
    """

    # Bad hack to circumvent t_coffee bug
    # Issue described already in: https://github.com/cbcrg/tcoffee/issues/3
    # Add an A in front of filename if the file begins with A
    filename=${test_alignment.fileName}
    first_letter_filename=\${filename:0:1}
    if [ "\$first_letter_filename" == "A" ]; then input="A"\$filename; else input=\$filename;  fi

    t_coffee -infile \$input -evaluate -output=score_ascii -lib \$lib -outfile "${test_alignment.baseName}.tcs"
    """
}

process SIM {
    container 'luisas/structural_regression:20'
    tag "SIM on $id - ${test_alignment.baseName}"
    storeDir "${params.outdir}/evaluation/sim/"
    label "process_low"

    input:
    tuple  val(id), file (test_alignment)

    output:
    path ("${test_alignment.baseName}.sim"), emit: sim

    script:
    """
    t_coffee -other_pg seq_reformat -in ${test_alignment} -output=sim > "${test_alignment.baseName}.sim"
    """
}

process EASEL_INFO {
    container 'edgano/hmmer:latest'
    tag "EASEL_INFO on $id"
    storeDir "${params.outdir}/evaluation/easel"

    input:
    tuple  val(id), file (test_alignment)


    output:
    tuple path("${test_alignment.baseName}.easel_INFO"), \
    path("${test_alignment.baseName}.avgLen"), \
    path("${test_alignment.baseName}.avgId"), emit: easelFiles

     shell:
     '''
     esl-alistat !{test_alignment} > !{test_alignment.baseName}.easel_INFO
     awk -F : '{ if (\$1=="Average length") printf "%s", \$2}' !{test_alignment.baseName}.easel_INFO | sed 's/ //g' > !{test_alignment.baseName}.avgLen
     awk -F : '{ if (\$1=="Average identity") printf "%s", substr(\$2, 1, length(\$2)-1)}' !{test_alignment.baseName}.easel_INFO | sed 's/ //g' > !{test_alignment.baseName}.avgId
     '''
}


process GAPS_PROGRESSIVE {
    container 'edgano/base:latest'
    tag "GAPS_PROG on $id"
    storeDir "${params.outdir}/evaluation/gaps"
    input:
    tuple  val (id), file (test_alignment)

    output:
    tuple val(id), \
    path("${test_alignment.baseName}.totGap"), \
    path("${test_alignment.baseName}.numSeq"), \
    path("${test_alignment.baseName}.alnLen"), emit: gapFiles

    script:
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    from decimal import *
    import os
    gap = '-'
    globalGap = 0
    avgGap = 0
    auxGap = 0
    totGapName= "${test_alignment.baseName}.totGap"
    numbSeqName= "${test_alignment.baseName}.numSeq"
    alnLenName= "${test_alignment.baseName}.alnLen"
    totGapFile= open(totGapName,"w+")
    numSeqFile= open(numbSeqName,"w+")
    alnLenFile= open(alnLenName,"w+")
    record = list(SeqIO.parse("${test_alignment}", "fasta"))
    for sequence in record:
        ## print(sequence.seq)
        auxGap = sequence.seq.count(gap)
        globalGap += auxGap
    avgGap = Decimal(globalGap) / Decimal(len(record))
    print "NumSeq: ",len(record)," GlobalGap: ",globalGap," AVG_Gap:",round(avgGap,3)
    totGapFile.write(str(globalGap))
    alnLenFile.write(str(len(record[0])))
    numSeqFile.write(str(len(record)))
    totGapFile.close()
    alnLenFile.close()
    numSeqFile.close()
    """
}


process EVAL_IRMSD{

  container 'luisas/structural_regression:7'
  tag "EVAL_ALIGNMENT on $id"
  storeDir "${params.outdir}/evaluation/irmsd/"
  label "process_low"

  input:
  tuple  val(id), file(msa), file(structures)

  output:
  path ("${msa.baseName}.total_irmsd"), emit: irmsd_summary

  script: 
  """
  # Prep templates
  for i in `awk 'sub(/^>/, "")' ${msa}`; do
      id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> template_list.txt
  done

  # Comp irmsd
  t_coffee -other_pg irmsd $msa -template_file template_list.txt | grep "TOTAL" > ${msa.baseName}.total_irmsd
  """
}

process EVAL_LIB{

  container 'luisas/structural_regression:7'
  tag "EVAL_ALIGNMENT_LIB on $id"
  storeDir  "${params.outdir}/evaluation/library/"
  label "process_low"

  input:
  tuple  val(id), file(msa), file(library)

  output:
  tuple val(id), path ("${msa.baseName}.score_ascii"), path("${msa.baseName}.html" ), emit: lib_summary

  script: 
  """
  filename=${msa}
  first_letter_filename=\${filename:0:1}
  if [ "\$first_letter_filename" == "A" ]; then input="A"${msa}; else input=${msa};  fi

  t_coffee -infile \$input  -lib $library -evaluate -output score_html -outfile ${msa.baseName}.html
  t_coffee -infile \$input  -lib $library -evaluate -output score_ascii -outfile ${msa.baseName}.score_ascii
  """
}




// iRMSD calculation


process EXTRACT_MSA_PAIRS {

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/${family_id}/alignment/pairs/"
  label 'process_small'
  tag "family_id"

  input:
  tuple val(family_id),file(msa)

  output:
  tuple val(family_id), file("${family_id}_pair*"), emit: pairs

  script:
  """
  python3 ${path_templates}/scripts/make_pairs.py $family_id $msa
  """
}

process GET_DISTANCES{

  container 'luisas/structural_regression:17'
  storeDir "${params.outdir}/${family_id}/evaluation/${evaluation_metric}/distances"
  label 'process_small'
  tag "$family_id"

  input:
  tuple val(family_id), file(msa), file(structures)
  each evaluation_metric
  each local_radius

  output:
  tuple val(family_id), file(msa), file("${msa.baseName}_${evaluation_metric}_${local_radius}.distances"), val(evaluation_metric), val(local_radius), emit: distances

  script:
  template "${path_templates}/scores/distances_${evaluation_metric}.sh"
}


process GET_SCORE{

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/${family_id}/evaluation/${evaluation_metric}"
  label 'process_small'
  tag "$family_id"

  input:
  tuple val(family_id), file(msa), file(distances), val(evaluation_metric), val(local_radius)

  output:
  tuple val(family_id),val(evaluation_metric), val(local_radius), file("${distances.baseName}_extended.score"), emit: score

  script:
  template "${path_templates}/scores/${evaluation_metric}.sh"
}

process MERGE_SCORES {

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/${family_id}/evaluation_summarized/"
  label 'process_small'
  tag "${family_id}_${evaluation_metric}_${local_radius}"

  input:
  tuple val(family_id), val(evaluation_metric), val(local_radius), file(score_files)

  output:
  tuple val(family_id), file("${family_id}_${evaluation_metric}_${local_radius}_scores.txt"), emit: scores_merged

  script:
  """
  for file in $score_files; do cat \$file >>  "${family_id}_${evaluation_metric}_${local_radius}_scores_temp.txt"; done

  # remove headers in the middle
  echo "aa1,aa2,id,irmsd,local_radius,res,tot_local_residues" > tmpfile
  sed -i '/^aa/d'  "${family_id}_${evaluation_metric}_${local_radius}_scores_temp.txt"
  cat "${family_id}_${evaluation_metric}_${local_radius}_scores_temp.txt" >> tmpfile
  mv tmpfile "${family_id}_${evaluation_metric}_${local_radius}_scores.txt"
  """
}
