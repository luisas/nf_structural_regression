
include { split_if_contains } from '../modules/functions.nf'

include { STRUCTURE_TO_3DI;   ENCODE_FASTA } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI } from '../modules/align.nf'
include { MERGE_MAPPINGS } from '../modules/utils.nf'
include { GENERATE_SEQUENCE_LIBRARY;} from '../modules/library.nf'
include { LIBRARY_3DI; GENERATE_3D_LIBRARY} from '../modules/library.nf'

workflow FOLDSEEK_LIBRARY {

  take:
    fastas
    structures
    matrix

  main:
    // Prepare the mapping file  
    fastas = fastas.map{ it -> [split_if_contains(it[0], "-ref", 0), it[0], it[1]]}
    STRUCTURE_TO_3DI(structures.groupTuple(by:0).map{it -> [it[0], it[2]]})
    MERGE_MAPPINGS(STRUCTURE_TO_3DI.out.mapping)
    mapping = MERGE_MAPPINGS.out.mapping
    fastas_to_map = fastas.combine(mapping, by:0)
    fastas_to_map = fastas_to_map.map{ it -> [it[1], it[2], it[3]]}
    // Prepare the fasta file
    ENCODE_FASTA(fastas_to_map)
    // Generate library 
    foldseek_library = LIBRARY_3DI(ENCODE_FASTA.out.encoded_fasta, matrix).library

  emit: 
    library = foldseek_library

}

workflow SEQUENCE_LIBRARY{

    take:
        fastas
    main:
    
        GENERATE_SEQUENCE_LIBRARY(fastas)
    
    emit: 
      library = GENERATE_SEQUENCE_LIBRARY.out.library
    
}


workflow TMALIGN_LIBRARY{

    take:
        seqs_and_trees_and_structures

    main:
        GENERATE_3D_LIBRARY(seqs_and_trees_and_structures, "TMalign_pair")
    
    emit: 
      library = GENERATE_3D_LIBRARY.out.library
    
}


workflow SAP_LIBRARY{

    take:
        seqs_and_trees_and_structures

    main:
        GENERATE_3D_LIBRARY(seqs_and_trees_and_structures, "sap_pair")
    
    emit: 
      library = GENERATE_3D_LIBRARY.out.library
    
}
