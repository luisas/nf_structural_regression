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
// target_db = Channel.fromPath( "${params.dbdir}/${params.target_db}",checkIfExists: true ).map { item -> [ item.baseName, item] }
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
         Target DB                                      : ${params.targetDB}
         """
         .stripIndent()

// import analysis pipelines
include { TREE_GENERATION } from './modules/treeGeneration'   params(params)
include { DYNAMIC_ANALYSIS } from './workflows/dynamic_analysis'    params(params)
include { REG_ANALYSIS } from './workflows/reg_analysis'        params(params)
include { PROG_ANALYSIS } from './workflows/prog_analysis'        params(params)
include { STRUCTURAL_REG_ANALYSIS } from './workflows/structural_reg_analysis'        params(params)
include { LIBRARIES_ANALYSIS } from './workflows/libraries_analysis'        params(params)
include { COMPACT_ANALYSIS } from './workflows/compact_analysis'        params(params)
include { MSA_3DI_ANALYSIS } from './workflows/msa_3di_analysis'        params(params)
include { REGFS_ANALYSIS } from './workflows/regfs_analysis'        params(params)
include { set_templates_path } from './modules/functions.nf'
path_templates = set_templates_path()

// Collect the FASTA files
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

// Collect the REFERENCE files
if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}

// Collect the STRUCTURES
print(params.structures_path)
str_path = params.structures_path
if(params.targetDB == "AF2_PRED"){
  structures_ch = Channel.fromPath(str_path)
                        .map { item -> [split_if_contains(item.getParent().getParent().baseName, "-ref", 0) , item.baseName.replace("_alphafold", ""), item] }
}else if(params.targetDB == "UniProtKB"){
  print("UniProtKB")
  structures_ch = Channel.fromPath(str_path)
                          .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0) , item.baseName.replace("_alphafold", ""), item] }
}
// if (params.alphafold) {
//   //str_path = "${params.af2_db_path}/colabfold/*/*_alphafold.pdb"
//   str_path = params.structures_path
//   // here was with af2 path, deep cleaning needed
//  // structures_ch = Channel.fromPath(str_path)
//  //                        .map { item -> [split_if_contains(item.getParent().getParent().baseName, "-ref", 0) , item.baseName.replace("_alphafold", ""), item] }

//   structures_ch = Channel.fromPath(str_path)
//                           .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0) , item.baseName.replace("_alphafold", ""), item] }


// }else{

//   structures_ch = Channel.fromPath("${params.experimental_structures_path}")
//                          .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0) , item.baseName.replace("_header", ""), item] }
// }


// Tokenize params
tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
library_method = params.library_methods.tokenize(',')
bucket_list = params.buckets.toString().tokenize(',')
dynamicX = params.dynamicX.toString().tokenize(',')
dynamicMasterAln = params.dynamicMasterAln.tokenize(',')
dynamicSlaveAln = params.dynamicSlaveAln.tokenize(',')
matrix = Channel.fromPath(params.matrix)
methodfile_path = "${path_templates}/method_files/${params.library_methods}.txt"
methodfile = Channel.fromPath(methodfile_path)


print(params.dataset_dir)
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
      PROG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method,structures_ch, target_db, channeled_dbs.collect())
    }

    if (params.dynamic_align){
      DYNAMIC_ANALYSIS(seqs_and_trees, refs_ch, tree_method, bucket_list, dynamicMasterAln, dynamicSlaveAln, dynamicX, structures_ch)
    }

    if (params.structural_regressive_align){
      STRUCTURAL_REG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list, target_db)
    }

    if (params.libraries_test){
      LIBRARIES_ANALYSIS(seqs_and_trees, refs_ch, library_method, tree_method, structures_ch, matrix)
    }
    
    if (params.compact_analysis){
      COMPACT_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }
    
    if (params.msafold){
      MSA_3DI_ANALYSIS(seqs_and_trees, refs_ch, library_method, tree_method, bucket_list, structures_ch, matrix)
    }

    if (params.regfold){
      print("Running the regressive fs analysis")
      REGFS_ANALYSIS(seqs_and_trees, refs_ch, library_method, tree_method, bucket_list, structures_ch, matrix, methodfile)
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
