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
 * Main XXX pipeline script
 *
 * @authors
 * Luisa Santus
 * Athanasios Baltzis
 * Edgar Garriga
 * Jose Espinosa-Carrasco
 */

//  example         https://github.com/nextflow-io/rnaseq-nf/tree/modules
/*
 * enables modules
 */
nextflow.enable.dsl = 2

/*
 * defaults parameter definitions
 */

//  Subsets of families - for testing
seq2improve="cryst,blmb,rrm,subt,ghf5,sdr,tRNA-synt_2b,zf-CCHH,egf,Acetyltransf,ghf13,p450,Rhodanese,aat,az,cytb,proteasome,GEL"
top20fam="gluts,myb_DNA-binding,tRNA-synt_2b,biotin_lipoyl,hom,ghf13,aldosered,hla,Rhodanese,PDZ,blmb,rhv,p450,adh,aat,rrm,Acetyltransf,sdr,zf-CCHH,rvp"
smallfam="seatoxin,scorptoxin,cyt3,rnasemam,bowman,hip,toxin,ghf11,TNF,sti"
testfam="seatoxin,bowman"

params.dataset_dir="/users/cn/lsantus/"
//params.dataset_dir="/home/luisasantus/Desktop/crg_cluster"
params.dataset = "homfam"
//params.dataset = "extHomfam_v35-uniprot"
params.seqs ="${params.dataset_dir}/data/structural_regression/${params.dataset}/combinedSeqs/{${testfam}}.fa"
params.refs = "${params.dataset_dir}/data/structural_regression/${params.dataset}/refs/{${testfam}}.ref"

//params.seqs ="${params.dataset_dir}/data/structural_regression/${params.dataset}/combinedSeqs/*.fa"
//params.refs = "${params.dataset_dir}/data/structural_regression/${params.dataset}/refs/*.ref"


params.align_methods = "FAMSA"
params.tree_methods = "FAMSA-medoid"

params.buckets = "100"
//  ## DYNAMIC parameters
params.dynamicX = "100000"
params.dynamicMasterAln="3dcoffee_msa"
params.dynamicMasterSize="100"
params.dynamicSlaveAln="famsa_msa"
params.dynamicSlaveSize="100000000"
params.dynamicConfig=true

if(params.dynamicMasterAln=="tcoffee_msa"){
  params.predict = true
}else{
  params.predict = true
}

params.cpu_flag="--cpu"
params.dynamic_align=true
params.regressive_align=true
params.progressive_align=false

params.evaluate=true
params.homoplasy=false
params.easel=false
params.metrics=false
params.compressAZ=false

// output directory
params.outdir = "$baseDir/results/${params.dataset}"
params.structures_dir="$baseDir/results/structures/"


log.info """\
         PIPELINE  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA))              : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
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
                  Dynamic config file                   : ${params.dynamicConfig}
                          master align - boundary       : ${params.dynamicMasterAln} - ${params.dynamicMasterSize}
                          slave align  - boundary       : ${params.dynamicSlaveAln} - ${params.dynamicSlaveSize}
                  Use AF2 predictions                   : ${params.predict}
         --##--
         Evaluate                                       : ${params.evaluate}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()

// import analysis pipelines
include { TREE_GENERATION } from './modules/treeGeneration'   params(params)
include { DYNAMIC_ANALYSIS } from './modules/dynamic_analysis'    params(params)
include { REG_ANALYSIS } from './modules/reg_analysis'        params(params)
include { PROG_ANALYSIS } from './modules/prog_analysis'        params(params)

// Channels containing sequences
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}

// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  trees = Channel.fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
}else {
  Channel.empty().set { trees }
}

// Tokenize params
tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
bucket_list = params.buckets.toString().tokenize(',')
dynamicX = params.dynamicX.toString().tokenize(',')

params.trees ="${params.dataset_dir}/data/structural_regression/${params.dataset}/trees/*/${tree_method}/{${testfam}}*.dnd"
//params.trees ="${params.dataset_dir}/data/structural_regression/${params.dataset}/trees/*/${tree_method}/*.dnd"

/*
 * main script flow
 */
workflow pipeline {

    TREE_GENERATION (seqs_ch, tree_method)
    seqs_ch
        .cross(TREE_GENERATION.out.trees)
        .map { it -> [ it[1][0], it[1][1], it[0][1], it[1][2] ] }
        .set { seqs_and_trees }

    // Run MSA
    if (params.regressive_align){
      REG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method, bucket_list)
    }

    if (params.progressive_align){
      PROG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method)
    }

    if (params.dynamic_align){
      DYNAMIC_ANALYSIS(seqs_and_trees, refs_ch, tree_method, bucket_list, dynamicX)
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

     sendMail(to: 'luisa.santus95@gmail.com', subject: 'My pipeline execution', body: msg)
 }
