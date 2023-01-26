#!/bin/bash nextflow
params.outdir = 'results'

process EVAL_ALIGNMENT {
    container 'luisas/structural_regression:7'
    tag "EVAL_ALIGNMENT on $id"
    storeDir "${params.outdir}/evaluation/score/"
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
    publishDir "${params.outdir}/evaluation/gaps", mode: 'copy', overwrite: true

    input:
    val align_type
    tuple  val (id), file (test_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path("*.totGap"), \
    path("*.numSeq"), \
    path("*.alnLen"), emit: gapFiles

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
    totGapName= "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.totGap"
    numbSeqName= "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.numSeq"
    alnLenName= "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.alnLen"
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


process TCS_optimize{

  container 'luisas/structural_regression:20'
  tag "TCS on $id - ${test_alignment.baseName}"
  storeDir "${params.outdir}/evaluation/tcs_optimization/"
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
