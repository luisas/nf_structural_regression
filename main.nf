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

/*
 * enables modules
 */
nextflow.enable.dsl = 2

/*
 * defaults parameter definitions
 */




//testfam="blm,egf,gpdh,lyase_1,int,subt,ldh,HLH,LIM,cyclo,proteasome,icd,msb,OTCace,HMG_box"
//testfam="seatoxin,scorptoxin"

testfam="seatoxin,scorptoxin"
//testfam="rub,ghf10,tgfb,sodcu,KAS,DMRL_synthase,tms,GEL,aadh,ace,ACPS,ARM,cat3,ccH,CH"
//testfam="kringle,cryst,DEATH,cah,mmp"
//testfam="Acetyltransf,sdr,rvp,zf-CCHH,rrm,aat,adh,p450,rhv,blmb,PDZ,Rhodanese"
testfammedium = "kringle,cryst,DEATH,cah,mmp,rub,ghf10,tgfb,sodcu,KAS,DMRL_synthase,tms,GEL"
testsmallhomfam="blm,egf,gpdh,lyase_1,int,subt,ldh,HLH,LIM,cyclo,proteasome,icd,msb,OTCace,HMG_box,flav,uce,peroxidase,sodfe,ghf1,cys,ace,glob,tim,hr,hormone_rec,hpr,oxidored_q6,asp,cytb,serpin,annexin,aadh,phc,ghf5,Ald_Xan_dh_2,mofe,Sulfotransfer,kunitz,GEL,tms,DMRL_synthase,KAS,sodcu,tgfb,ghf10,rub,mmp,cah,DEATH,cryst,kringle,az,il8,ltn,phoslip,slectin,trfl,ins,ChtBD,ghf22,ricin,profilin,Stap_Strp_toxin,sti,TNF,ghf11,toxin,bowman,rnasemam,cyt3,scorptoxin,hip,seatoxin,test"
testfammedium2="rrm,aat,adh,p450,rhv,blmb,PDZ,Rhodanese,hla,aldosered,ghf13,hom,biotin_lipoyl,tRNA-synt_2b,myb_DNA-binding,gluts"
testlarge2="Acetyltransf,sdr,rvp,zf-CCHH"

testfamlarge= "aadh,ace,ACPS,ARM,cat3,ccH,CH"
testfamxlarge= "aabp,actin,adk,C2,cox,COX2,CPS"
testfamhuge= "AAA,ABC_tran,blmb,fn3,rep,sdr"

params.dataset_dir="/users/cn/lsantus/"
params.dataset = "homfam"
fasta_dirs="combinedSeqs,refs"
params.seqs ="${params.dataset_dir}/data/structural_regression/${params.dataset}/{$fasta_dirs}/*.fa"
params.refs = "${params.dataset_dir}/data/structural_regression/${params.dataset}/refs/*.ref"
params.af2_db_path = "${params.dataset_dir}/data/structural_regression/af2_structures"

params.align_methods = "FAMSA"
params.tree_methods = "FAMSA-medoid"

params.buckets = "50,1000"
//  ## DYNAMIC parameters
params.dynamicX = "1"
params.dynamicMasterAln="famsa_msa"
params.dynamicSlaveAln="famsa_msa"

params.predict = true


params.dynamic_align=false
params.regressive_align=true
params.progressive_align=true

params.evaluate=true


// output directory
params.outdir = "$baseDir/results/${params.dataset}"


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
         """
         .stripIndent()

// import analysis pipelines
include { TREE_GENERATION } from './modules/treeGeneration'   params(params)
include { DYNAMIC_ANALYSIS } from './modules/dynamic_analysis'    params(params)
include { REG_ANALYSIS } from './modules/reg_analysis'        params(params)
include { PROG_ANALYSIS } from './modules/prog_analysis'        params(params)

// Channels
seqs_ch = Channel.fromPath( params.seqs, checkIfExists: true ).map { item -> [ item.baseName, item] }

seqs_ch.view()


structures_ch = Channel.fromPath("${params.af2_db_path}/colabfold_header/*/**/*.pdb")
                       .map { item -> [ item.getParent().getParent().baseName, item.baseName, item] }

if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}


// Tokenize params
tree_method = params.tree_methods.tokenize(',')
align_method = params.align_methods.tokenize(',')
bucket_list = params.buckets.toString().tokenize(',')
dynamicX = params.dynamicX.toString().tokenize(',')
dynamicMasterAln = params.dynamicMasterAln.tokenize(',')
dynamicSlaveAln = params.dynamicSlaveAln.tokenize(',')


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
      PROG_ANALYSIS(seqs_and_trees, refs_ch, align_method, tree_method)
    }

    if (params.dynamic_align){
      DYNAMIC_ANALYSIS(seqs_and_trees, refs_ch, tree_method, bucket_list, dynamicMasterAln, dynamicSlaveAln, dynamicX, structures_ch)
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
