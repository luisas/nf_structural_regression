#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.dataset_dir="/users/cn/lsantus"
include { split_if_contains } from './modules/functions.nf'
include { MMSEQS_SEARCH;TEMPLATE_FROM_DB_HITS;FETCH_STRUCTURES; GET_GDT; ADD_HEADER; PREP_STRUCTURES; FETCH_FASTA; RENAME_PDB; FETCH_STRUCTURES_UNIPROT }        from './modules/structures.nf'


//params.dataset_dir="/home/luisasantus/Desktop/crg_cluster"
params.dataset = "homfam"
//params.fastas = "${params.dataset_dir}/data/structural_regression/${params.dataset}/refs/*.fa"
params.fastas = "${params.dataset_dir}/data/structural_regression/${params.dataset}/refs/{$params.testfam}*.fa"


//params.target_db = "PDB"
params.target_db = "UniProtKB"

params.dbdir = "${params.dataset_dir}/data/db"
//params.outdir = "${params.dataset_dir}/data/structural_regression/${params.dataset}/pdbs"


// Params for search
target_db = Channel.fromPath( "${params.dbdir}/${params.target_db}",checkIfExists: true ).map { item -> [ item.baseName, item] }
// Here i use the nulls so that the cardinality of the Channel matches the one of the previous
refs_ch = Channel.fromPath( "$params.fastas" ).map { item -> [ item.baseName,null,item,null] }

println("$params.target_db")

//af2_structures = Channel.fromPath( "$params.af2_structures" ).map{ it -> [ it.getParent().getParent().baseName, it.baseName, it ] }
//af2_structures.view()


workflow prep_structures_pdb {

  // 1. Find sequence hits in PDB
  MMSEQS_SEARCH(refs_ch,target_db.collect())
  // 2. Create the template file and obtain the best mmseqs hit
  // The mmseqs hit comes with start and end positions of the match, which I will later cut.
  mmseqs_hits = TEMPLATE_FROM_DB_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 })
  // 3. Download the structures
  //    and cut them according to the positions the hits (extract only the real matching chunk)
  FETCH_STRUCTURES(mmseqs_hits)
  fetched_structures = FETCH_STRUCTURES.out.fetched_structures.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[2], it[3], it[4]] }
                                         .transpose()
                                         .map{ it -> [ it[0], split_if_contains(it[3].baseName, "_ref", 0), it[1], it[2], it[3] ]}
                                        // .filter { it[1] =~/.*_ref/ }.map{ it -> [ it[0], split_if_contains(it[1], "_ref", 0), it[2]]}

  // 3b. Download the fastas
  FETCH_FASTA(mmseqs_hits)
  fetched_fastas = FETCH_FASTA.out.fastas.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[2], it[3], it[4]] }
                                         .transpose()
                                         .map{ it -> [ it[0], split_if_contains(it[3].baseName, "_ref", 0), it[1], it[2], it[3] ]}

  structures = fetched_structures.join(fetched_fastas, by: [0,1]).map{ it -> [ it[0], it[1], it[2], it[3], it[4], it[7]]}

  // 4. Prep PDB (extract_from_pdb)
  PREP_STRUCTURES(structures)
  //PREP_STRUCTURES.out.structures.view()
  preprocessed_structures = PREP_STRUCTURES.out.structures.map{ it -> [ it[0], it[1], it[4]]}
  //preprocessed_structures.view()
  reference_structures = ADD_HEADER(preprocessed_structures)
  RENAME_PDB(reference_structures)

}

workflow prep_structures_uniprot {

  // 1. Find sequence hits in PDB
  MMSEQS_SEARCH(refs_ch,target_db.collect())
  // 2. Create the template file and obtain the best mmseqs hit
  // The mmseqs hit comes with start and end positions of the match, which I will later cut.
  mmseqs_hits = TEMPLATE_FROM_DB_HITS(MMSEQS_SEARCH.out.hits.filter{ it[2].size()>0 })
  // 3. Download the structures
  //    and cut them according to the positions the hits (extract only the real matching chunk)
  FETCH_STRUCTURES_UNIPROT(mmseqs_hits)
  fetched_structures = FETCH_STRUCTURES_UNIPROT.out.fetched_structures.map{ it -> [split_if_contains(it[0], "-ref", 0) , it[2], it[3], it[4]] }
                                         .transpose()
                                         .map{ it -> [ it[0], split_if_contains(it[3].baseName, "_ref", 0), it[1], it[2], it[3] ]}
                                        // .filter { it[1] =~/.*_ref/ }.map{ it -> [ it[0], split_if_contains(it[1], "_ref", 0), it[2]]}
}


workflow {
  //prep_structures_pdb()
  prep_structures_uniprot()
}
