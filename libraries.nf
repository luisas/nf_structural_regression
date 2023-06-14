#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { split_if_contains } from './modules/functions.nf'
include { TREE_GENERATION } from './modules/treeGeneration'   params(params)
include { MODIFY_LIBRARY; MERGE_LIBRARIES; AGGREGATE_LIBRARY; MERGE_LIBRARIES_MANUAL} from './modules/library'   params(params)
include { set_templates_path } from './modules/functions.nf'
path_templates = set_templates_path()
include { EVALUATE_MSA ; EVALUATE_MSA_STRUCTURAL; EVALUATE_MSA_LIBRARIES;  } from './subworkflows/evaluate.nf'
include { EVALUATE_MSA as EVALUATE_MSA_2 ; EVALUATE_MSA_STRUCTURAL as EVALUATE_MSA_STRUCTURAL_2 ; EVALUATE_MSA_LIBRARIES as EVALUATE_MSA_LIBRARIES_2;  } from './subworkflows/evaluate.nf'
include { EVALUATE_MSA as EVALUATE_MSA_3 ; EVALUATE_MSA_STRUCTURAL as EVALUATE_MSA_STRUCTURAL_3 ; EVALUATE_MSA_LIBRARIES as EVALUATE_MSA_LIBRARIES_3;  } from './subworkflows/evaluate.nf'
include { EVALUATE_MSA as EVALUATE_MSA_4 ; EVALUATE_MSA_STRUCTURAL as EVALUATE_MSA_STRUCTURAL_4 ; EVALUATE_MSA_LIBRARIES as EVALUATE_MSA_LIBRARIES_4;  } from './subworkflows/evaluate.nf'

process ALIGN_WITH_LIBRARY{
    container 'luisas/structural_regression:20'
    storeDir "${params.outdir}/alignment/$id/"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(lib_method), path(library)


    output:
    tuple val(id), path("${library.baseName}.aln"), emit: alignment

    script:
    """
    t_coffee -infile $seqs -lib $library -usetree=$guide_tree -outfile ${library.baseName}.aln -output fasta_aln
    """

}

process ALIGN_WITH_LIBRARY_2{
    container 'luisas/structural_regression:20'
    storeDir "${params.outdir}/alignment/$id/"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(lib_method), path(library)


    output:
    tuple val(id), path("${library.baseName}.aln"), emit: alignment

    script:
    """
    t_coffee -infile $seqs -lib $library -usetree=$guide_tree -outfile ${library.baseName}.aln -output fasta_aln
    """

}

process ALIGN_WITH_LIBRARY_3{
    container 'luisas/structural_regression:20'
    storeDir "${params.outdir}/alignment/$id/"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(lib_method), path(library)


    output:
    tuple val(id), path("${library.baseName}.aln"), emit: alignment

    script:
    """
    t_coffee -infile $seqs -lib $library -usetree=$guide_tree -outfile ${library.baseName}.aln -output fasta_aln
    """

}


process ALIGN_WITH_LIBRARY_4{
    container 'luisas/structural_regression:20'
    storeDir "${params.outdir}/alignment/$id/"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), val(lib_method), path(library)


    output:
    tuple val(id), path("${library.baseName}.aln"), emit: alignment

    script:
    """
    t_coffee -infile $seqs -lib $library -usetree=$guide_tree -outfile ${library.baseName}.aln -output fasta_aln
    """

}

seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }
libaries = Channel.fromPath( params.libs, checkIfExists: true ).map { item -> [ item.getParent().baseName, item.getParent().getParent().getParent().baseName, item] }
seq_lib_ch = Channel.fromPath( params.seq_lib, checkIfExists: true ).map { item -> [ item.getParent().baseName, item] }

tree_method = params.tree_methods.tokenize(',')

mod_methods = params.mod_methods.tokenize(',')
N = params.Ns.tokenize(',')
refs = Channel.fromPath( params.refs, checkIfExists: true ).map { item -> [ item.baseName, item] }
structures = Channel.fromPath(params.structures_path)
                        .map { item -> [split_if_contains(item.getParent().getParent().baseName, "-ref", 0) , item.baseName.replace("_alphafold", ""), item] }

workflow pipeline {

 TREE_GENERATION (seqs_ch, tree_method)
    seqs_ch
        .cross(TREE_GENERATION.out.trees)
        .map { it -> [ it[1][0], it[1][1], it[0][1], it[1][2] ] }
        .set { seqs_and_trees }
    seqs_trees_libraries = seqs_and_trees.combine(libaries, by:0)


    // 1 . ALIGN WITH THE NORMAL LIBRARY 
    ALIGN_WITH_LIBRARY(seqs_trees_libraries)
    EVALUATE_MSA( ALIGN_WITH_LIBRARY.out.alignment, refs)
    EVALUATE_MSA_STRUCTURAL ( ALIGN_WITH_LIBRARY.out.alignment, structures)

    // 2. MODIFY THE LIBRARY AND ALIGN WITH THE MODIFIED LIBRARY
    if (params.modify){
        MODIFY_LIBRARY (seqs_trees_libraries, mod_methods, N)
        ALIGN_WITH_LIBRARY_2(MODIFY_LIBRARY.out.lib)
        EVALUATE_MSA_2( ALIGN_WITH_LIBRARY_2.out.alignment, refs)
        EVALUATE_MSA_STRUCTURAL_2 ( ALIGN_WITH_LIBRARY_2.out.alignment, structures)
    }

    // 3. MERGE THE MODIFIED LIBRARY AND THE SEQUENCE LIBRARY. ALIGN WITH THE MERGED LIBRARY
    if(params.merge_auto){
        scaled_lib = MODIFY_LIBRARY.out.lib.map{ it ->  [it[0], it[5]]}
        libs = scaled_lib.combine(seq_lib_ch, by:0)
        libs.map{ it -> [it[0]+"_sequence_${params.lib}_${params.mod_methods}_${params.Ns}", it[1], it[2]]}.set{libs}
        MERGE_LIBRARIES(libs, "sequence_${params.lib}_${params.mod_methods}_${params.Ns}")
        MERGE_LIBRARIES.out.library.map{ it ->  [it[0].split("_")[0], it[0]+"_${params.aggfunc}", it[1]]}.set{merged_lib}
        print("${params.aggfunc}")
        AGGREGATE_LIBRARY(merged_lib, "${params.aggfunc}")
        ALIGN_WITH_LIBRARY_3(seqs_and_trees.combine(AGGREGATE_LIBRARY.out.lib, by:0))
        ALIGN_WITH_LIBRARY_3.out.alignment.view()
        EVALUATE_MSA_3( ALIGN_WITH_LIBRARY_3.out.alignment, refs)
        EVALUATE_MSA_STRUCTURAL_3 ( ALIGN_WITH_LIBRARY_3.out.alignment, structures)
    }
    else if(params.merge_manual){
        scaled_lib = MODIFY_LIBRARY.out.lib.map{ it ->  [it[0], it[5]]}
        libs = scaled_lib.combine(seq_lib_ch, by:0)
        libs.map{ it -> [it[0]+"_sequence_${params.lib}_${params.mod_methods}_${params.Ns}", it[1], it[2]]}.set{libs}
        MERGE_LIBRARIES_MANUAL(libs, "sequence_${params.lib}_${params.mod_methods}_${params.Ns}", "${params.aggfunc}")
        MERGE_LIBRARIES_MANUAL.out.library.map{ it ->  [it[0].split("_")[0], it[0]+"_${params.aggfunc}", it[2]]}.set{merged_lib}
        merged_lib.view()
        seqs_and_trees.combine(merged_lib, by:0).view()
        ALIGN_WITH_LIBRARY_4(seqs_and_trees.combine(merged_lib, by:0))
        ALIGN_WITH_LIBRARY_4.out.alignment.view()
        EVALUATE_MSA_4( ALIGN_WITH_LIBRARY_4.out.alignment, refs)
        EVALUATE_MSA_STRUCTURAL_4 ( ALIGN_WITH_LIBRARY_4.out.alignment, structures)
    }


}
workflow {
  pipeline()
}