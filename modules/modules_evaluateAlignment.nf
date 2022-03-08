#!/bin/bash nextflow
params.outdir = 'results'

process EVAL_COMPRESS {
    container 'edgano/tcoffee:pdb'
    tag "EVAL_COMPRESS on $id"
    publishDir "${params.outdir}/score/tc", mode: 'copy', overwrite: true, pattern: '*.tc'
    publishDir "${params.outdir}/score/sp", mode: 'copy', overwrite: true, pattern: '*.sp'
    publishDir "${params.outdir}/score/col", mode: 'copy', overwrite: true, pattern: '*.col'

    input:
    val align_type
    tuple  val (id), file (test_alignment), file (ref_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"), emit: spScore

    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"), emit: tcScore

    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"), emit: colScore

    script:
    """
    ## extract ids
    grep \"^>\"  ${ref_alignment} | sed 's/^>//g' > id_extracted.txt
    ## extract embedded ref from aln
    awk 'NR==FNR{n[">"\$0];next} f{print f ORS \$0;f=""} \$0 in n{f=\$0}' id_extracted.txt ${test_alignment} > ref_extracted.aln
    ##  decompress ref_aln, from gap format to fasta
    ## TODO -> place it in the container.
    gcc ${baseDir}/bin/bioinfoCommands/gap.c -o gap
    chmod 777 gap
    ./gap ref_extracted.aln false                ## resultG2F.fa
    ##  compare ref and ref_aln
        ## Sum-of-Pairs Score ##
        t_coffee -other_pg aln_compare \
                -al1 ${ref_alignment} \
                -al2 resultG2F.fa \
                -compare_mode sp \
                | grep -v "seq1" | grep -v '*' | \
                awk '{ print \$4}' ORS="\t" \
                > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"
        ## Total Column Score ##
        t_coffee -other_pg aln_compare \
                -al1 ${ref_alignment} \
                -al2 resultG2F.fa \
                -compare_mode tc \
                | grep -v "seq1" | grep -v '*' | \
                awk '{ print \$4}' ORS="\t" \
                > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"
        ## Column Score ##
        t_coffee -other_pg aln_compare \
                -al1 ${ref_alignment} \
                -al2 resultG2F.fa \
                -compare_mode column \
                | grep -v "seq1" | grep -v '*' | \
                awk '{ print \$4}' ORS="\t" \
                > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"
    """
}

process EVAL_ALIGNMENT {
    container 'edgano/tcoffee:pdb'
    tag "EVAL_ALIGNMENT on $id"
    publishDir "${params.outdir}/score/tc", mode: 'copy', overwrite: true, pattern: '*.tc'
    publishDir "${params.outdir}/score/sp", mode: 'copy', overwrite: true, pattern: '*.sp'
    publishDir "${params.outdir}/score/col", mode: 'copy', overwrite: true, pattern: '*.col'

    input:
    val align_type
    tuple  val (id), file (test_alignment), file (ref_alignment)
    val align_method
    val tree_method
    val bucket_size

    output:
    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"), emit: spScore

    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"), emit: tcScore

    tuple val(id), \
    val(align_type), \
    val(bucket_size), \
    val(align_method), \
    val(tree_method), \
    path ("${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"), emit: colScore

    script:
    """
    ## Sum-of-Pairs Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.sp"
    ## Total Column Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.tc"
    ## Column Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode column \
            | grep -v "seq1" | grep -v '*' | \
              awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.with.${tree_method}.tree.col"
    """
}
