#!/bin/bash nextflow
include { FOLDSEEK_LIBRARY; SEQUENCE_LIBRARY; TMALIGN_LIBRARY;SAP_LIBRARY } from '../subworkflows/prep_libraries.nf'
include { split_if_contains } from '../modules/functions.nf'
include { COMPACT_ALIGNER; FS_ALIGNER; FSREG_ALIGNER} from '../modules/align.nf'
include { STRUCTURE_TO_3DI;   ENCODE_FASTA } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI } from '../modules/align.nf'
include { MERGE_MAPPINGS } from '../modules/utils.nf'
include { PREP_FS_SEQS } from '../modules/encoding.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'

workflow REGFS_ANALYSIS {

  take:
    seqs_and_trees
    refs
    library_method
    tree_method
    bucket_size
    structures
    matrix
    methodfile

  main:

    // Prep channels 
    seqs_and_trees = seqs_and_trees.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1], it[2], it[3]]}
    seqs_and_trees_and_structures = seqs_and_trees.combine(structures, by: [0]).groupTuple(by:[1,2,3,4])
                                                    .map { it -> [ it[1], it[2], it[3], it[4], it[6]]}
    fastas = seqs_and_trees.map { it -> [it[1], it[3]]}

    // Prep sequences foldseek
    fastas = fastas.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1]]}
    STRUCTURE_TO_3DI(structures.groupTuple(by:0).map{it -> [it[0], it[2]]})
    MERGE_MAPPINGS(STRUCTURE_TO_3DI.out.mapping)
    mapping = MERGE_MAPPINGS.out.mapping
    fastas_to_map = fastas.combine(mapping, by:0)
    fastas_to_map = fastas_to_map.map{ it -> [it[1], it[2], it[3]]}
    PREP_FS_SEQS(fastas_to_map)
    seqs_and_trees_aln = seqs_and_trees.map{ it -> [it[1], it[2], it[3], it[4]]}
    seqs_and_trees_aln2 = seqs_and_trees_aln.combine(PREP_FS_SEQS.out.fsdir, by:0)

    FSREG_ALIGNER(seqs_and_trees_aln2, matrix.collect(), library_method, methodfile.collect(), bucket_size )
    if (params.evaluate){
             EVALUATE_MSA( FSREG_ALIGNER.out.alignmentFile, refs)
    }

}

