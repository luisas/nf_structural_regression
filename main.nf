#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Main Structural regression  pipeline script
 *
 * @authors
 * Luisa Santus
 * Athanasios Baltzis
 * Edgar Garriga
 * Jose Espinosa-Carrasco
 */

/*
 * enables modules
 */
nextflow.enable.dsl = 2

/*
 * defaults parameter definitions
 */



include { split_if_contains } from './modules/functions.nf'

// Target database (MMSEQS)
target_db = Channel.fromPath( "${params.dbdir}/${params.target_db}",checkIfExists: true ).map { item -> [ item.baseName, item] }
// Target database (PDB)
channeled_dbs = Channel.fromPath( params.blast_database ).map { item -> [only_first_extension(item.name), item] }.groupTuple( by:[0]).map{ it -> [it[0], it[1][0],it[1][1..-1]]}

log.info """\
         STRUCTURAL REGRESSION  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA))              : ${params.refs}
         Input trees (NEWICK)                           : ${params.outdir}/trees/
         Alignment methods                              : ${params.align_methods}
         Tree methods                                   : ${params.tree_methods}
         Bucket size                                    : ${params.buckets}
         --##--
         Generate Progressive alignments                : ${params.progressive_align}
         --##--
         Generate Regressive alignments                 : ${params.regressive_align}
         --##--
         Generate Dynamic alignments                    : ${params.dynamic_align}
                  Dynamic size                          : ${params.dynamicX}
                  Dynamic config
                          master align                  : ${params.dynamicMasterAln}
                          slave align                   : ${params.dynamicSlaveAln}
                  Use AF2 predictions                   : ${params.predict}
         --##--
         Evaluate                                       : ${params.evaluate}
         Output directory (DIRECTORY)                   : ${params.outdir}
         --##--
         Target DB                                      : ${params.dbdir}
         Blast DB                                      : ${params.blast_database}
         """
         .stripIndent()

// import analysis pipelines
include { TREE_GENERATION } from './modules/treeGeneration'   params(params)
include { DYNAMIC_ANALYSIS } from './workflows/dynamic_analysis'    params(params)
include { REG_ANALYSIS } from './workflows/reg_analysis'        params(params)
include { PROG_ANALYSIS } from './workflows/prog_analysis'        params(params)
include { STRUCTURAL_REG_ANALYSIS } from './workflows/structural_reg_analysis'        params(params)
include { FOLDSEEK_ANALYSIS } from './workflows/foldseek_analysis'        params(params)
include { LIBRARIES_ANALYSIS } from './workflows/libraries_analysis'        params(params)
include { COMPACT_ANALYSIS } from './workflows/compact_analysis'        params(params)

// Collect the FASTA files
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

// Collect the STRUCTURES
if (params.alphafold) {
  //str_path = "${params.af2_db_path}/colabfold/*/*_alphafold.pdb"
  str_path = params.structures_path
  structures_ch = Channel.fromPath(str_path)
                         .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0) , item.baseName.replace("_alphafold", ""), item] }
}else{

  structures_ch = Channel.fromPath("${params.experimental_structures_path}")
                         .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0) , item.baseName.replace("_header", ""), item] }
}

if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}

//print("${params.af2_db_path}")
//structures_ch.view()

// Tokenize params
tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
if(params.libraries_test){
  library_method = params.libraries_test.tokenize(',')
}
bucket_list = params.buckets.toString().tokenize(',')
dynamicX = params.dynamicX.toString().tokenize(',')
dynamicMasterAln = params.dynamicMasterAln.tokenize(',')
dynamicSlaveAln = params.dynamicSlaveAln.tokenize(',')

if(params.libraries_test){
  matrix = Channel.fromPath(params.matrix)
}
/*
 * main script flow
 */
workflow pipeline {

    TREE_GENERATION (seqs_ch, tree_method)
    seqs_ch
        .cross(TREE_GENERATION.out.trees)
        .map { it -> [ it[1][0], it[1][1], it[0][1], it[1][2] ] }
        .set { seqs_and_trees }

    if (params.regressive_align){
      REG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }

    if (params.progressive_align){
      // here the situation about the structures need to be either removed or polished!
      PROG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method,structures_ch, target_db, channeled_dbs.collect())
    }

    if (params.dynamic_align){
      DYNAMIC_ANALYSIS(seqs_and_trees, refs_ch, tree_method, bucket_list, dynamicMasterAln, dynamicSlaveAln, dynamicX, structures_ch)
    }

    if (params.structural_regressive_align){
      STRUCTURAL_REG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list, target_db)
    }

    if (params.foldseek_align){
      FOLDSEEK_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, structures_ch)
    }

    if (params.libraries_test){
      LIBRARIES_ANALYSIS(seqs_and_trees, refs_ch, library_method, tree_method, structures_ch, matrix)
    }
    
    if (params.compact_analysis){
      COMPACT_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }
}

workflow {
  pipeline()
}

/*
 * completion handler
 */
 workflow.onComplete {

     def msg = """\
         Pipeline execution summary
         ---------------------------
         Completed at: ${workflow.complete}
         Duration    : ${workflow.duration}
         Success     : ${workflow.success}
         workDir     : ${workflow.workDir}
         exit status : ${workflow.exitStatus}
         """
         .stripIndent()

     sendMail(to: 'luisa.santus95@gmail.com', subject: 'Structural regression pipeline', body: msg)
 }
 /*
  * Extra functions
  */

 def only_first_extension(String filename){
   return filename.split("\\.").take(2).join('.')
 }
