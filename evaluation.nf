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


nextflow.enable.dsl = 2

// import analysis pipelines
include { split_if_contains } from './modules/functions.nf'
include { set_templates_path } from './modules/functions.nf'
path_templates = set_templates_path()
include { EVALUATE_MSA_STRUCTURAL } from './subworkflows/evaluate.nf'

// Collect the ALIGNMNTS files
alignments_ch = Channel.fromPath( params.alignments, checkIfExists: true ).map { item -> [ item.baseName.split("\\.")[0], item] }

// Collect the REFERENCE files
if ( params.refs ) {
  refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }
}

// Collect the STRUCTURES
if(params.targetDB == "UniProtKB"){
  structures_ch = Channel.fromPath(params.structures_path)
                          .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0), item] }
                          .groupTuple(by:0)
}else if(params.targetDB == "AF2_PRED"){
  structures_ch = Channel.fromPath(params.structures_path)
                          .map { it -> [split_if_contains(it.getParent().getParent().baseName, "-ref", 0), it] }
                          .groupTuple(by:0)
}

/*
 * main script flow
 */
workflow pipeline {

    if (params.structural){
      EVALUATE_MSA_STRUCTURAL(alignments_ch, structures_ch)
    }


}



workflow {
  pipeline()
}



