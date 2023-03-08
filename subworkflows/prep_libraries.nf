
include { split_if_contains } from '../modules/functions.nf'

include { STRUCTURE_TO_3DI;   ENCODE_FASTA } from '../modules/encoding.nf'
include { ALIGN_WITH_3DI } from '../modules/align.nf'
include {MERGE_MAPPINGS } from '../modules/utils.nf'
include {GENERATE_3DI_LIBRARY; GENERATE_SEQUENCE_LIBRARY} from '../modules/library.nf'

workflow FOLDSEEK_LIBRARY {

  take:
    fastas
    structures
    matrix

  main:

    // Prepare the mapping file  
    STRUCTURE_TO_3DI(structures.groupTuple(by:0).map{it -> [it[0], it[2]]})
    MERGE_MAPPINGS(STRUCTURE_TO_3DI.out.mapping)
    mapping = MERGE_MAPPINGS.out.mapping
    fastas_to_map = fastas.combine(mapping, by:0)

    // Prepare the fasta file
    ENCODE_FASTA(fastas_to_map)
    // Generate library 
    GENERATE_3DI_LIBRARY(ENCODE_FASTA.out.encoded_fasta, matrix)
    // TODO: change back the sequence to original
    foldseek_library = GENERATE_3DI_LIBRARY.out.library

  // emit: 
  //  library = foldseek_library

}

workflow SEQUENCE_LIBRARY{

    take:
        fastas
    main:
    
        GENERATE_SEQUENCE_LIBRARY(fastas)
    
    emit: 
      library = GENERATE_SEQUENCE_LIBRARY.out.library
    
}