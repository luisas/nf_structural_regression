#!/bin/bash nextflow
params.outdir = 'results'

include {GENERATE_DYNAMIC_CONFIG; EXTRACT_SEQUENCES; ADD_PDB_HEADERS; CHECK_CACHE}      from './preprocess.nf'
include {ADD_PDB_HEADERS as ADD_PDB_HEADERS_2}      from './preprocess.nf'

include {EVAL_ALIGNMENT}      from './modules_evaluateAlignment.nf'
include {PROG_ALIGNER; PROG_ALIGNER_STRUCTURES}       from './generateAlignment.nf'
include {EASEL_INFO}        from './modules_evaluateAlignment.nf'
include { RUN_COLABFOLD } from './localcolabfold.nf'
include { split_if_contains } from './functions.nf'
include {MMSEQS_PREP_DB; MMSEQS_SEARCH; TEMPLATE_FROM_DB_HITS; FETCH_STRUCTURES}        from './structures.nf'


workflow PROG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    structures
    target_db

  main:

    align_method_string = align_method.toString().replace("[", "").replace("]", "")
    println(align_method_string)

    if(align_method_string == "EXPRESSO"){

      // Collect the structures with
      MMSEQS_SEARCH(seqs_and_trees,target_db.collect())
      TEMPLATE_FROM_DB_HITS(MMSEQS_SEARCH.out.hits)
      FETCH_STRUCTURES(TEMPLATE_FROM_DB_HITS.out.filtered_hits)
      seqs_trees_structures = seqs_and_trees.combine(FETCH_STRUCTURES.out.fetched_structures, by: [0])

      PROG_ALIGNER_STRUCTURES (seqs_trees_structures, align_method)
      prog_alignment = PROG_ALIGNER_STRUCTURES.out.alignmentFile


    }else if(align_method_string == "3DCOFFEE"){

      /// HERE MAJOR CODE DUPLICATION WITH DYNAMIC - IF WE WANT TO KEEP, ORDER!
      ids_done = structures.ifEmpty(">--")
                           .collectFile() { item -> [ "ids_done.txt", item[1] + '\n' ]}
                           .collect()

      seq_and_trees_mod = seqs_and_trees.map{ it -> [it[0], it[1], "MAX", it[2]]}

      CHECK_CACHE(seq_and_trees_mod, ids_done)
      RUN_COLABFOLD(CHECK_CACHE.out.seqToPredict.filter{ it[3].size()>0 }.splitFasta( by: params.n_af2, file: true ),
                    params.model_type,
                    params.db )


      precomputed_structures = CHECK_CACHE.out.idsDone
                                          .map{ item -> [item.baseName.split("_")[0..-3].join("_"),
                                                         item.baseName.split("_")[-2],
                                                         item.baseName.split("_")[-1],
                                                         item.splitText().collect { it.trim()},
                                                         ]}
                                          .transpose()
                                          .map{ it -> [ split_if_contains(it[0], "-ref", 0), it[3].replace("/", "_"), it[1], it[2], it[0]]}
                                          .join(structures, by: [0,1])
                                          .map{ it -> [it[4], it[2], it[3], it[5]]}
      // 2. Get a channel with both the newly and pre- computed structures
      all_structures = RUN_COLABFOLD.out.af2_pdb.groupTuple(by:[0,1,2])
                                    .map{ it -> [it[0],it[1],it[2],it[3].flatten()]}
                                    .concat(precomputed_structures)
                                    .map{ it -> [it[0], it[1], it[3].toString().split("/")[-1].split("_alphafold")[0], it[3]] }


      // How it looks like: all structures -->  [fam,tree,sequence_id, sequence_id.pdb]
      ADD_PDB_HEADERS(all_structures.map{ it -> [ it[0], it[1], it[2], it[3], split_if_contains(it[0], "-ref", 0)]})

      // Combine by family and tree
      seqs_and_trees
        .combine(ADD_PDB_HEADERS.out.pdb.groupTuple(by: [0,1]), by: [0,1])
        .set{ seqs_and_trees_and_structures }


      PROG_ALIGNER_STRUCTURES (seqs_and_trees_and_structures, align_method)
      prog_alignment = PROG_ALIGNER_STRUCTURES.out.alignmentFile



    }else{

      PROG_ALIGNER (seqs_and_trees, align_method)
      prog_alignment = PROG_ALIGNER.out.alignmentFile

    }


    /*
    *   EVALUATION
    */
    if (params.evaluate){

      alignments = prog_alignment.map{ it -> [ it[0].replaceAll("-ref", ""), it[1] ] }

      refs_ch
        .cross (alignments)
        .map { it -> [ it[1][0], it[1][1], it[0][1] ] }.view()
        .set { alignment_and_ref }
      EVAL_ALIGNMENT (alignment_and_ref)
      EASEL_INFO (alignment_and_ref)

      EVAL_ALIGNMENT.out.scores
                    .collectFile(name: "progressive.scores.csv", newLine: true, storeDir:"${params.outdir}/evaluation/CSV/")
    }

    emit:
    alignment = prog_alignment
}
