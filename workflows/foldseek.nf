
include { STRUCTURE_TO_3DI; MERGE_MAPPINGS; PREP_FS_SEQS } from '../modules/encoding.nf'


workflow FOLDSEEK_CONVERT{
    take:
        structures

    main:
        STRUCTURE_TO_3DI(structures)
        MERGE_MAPPINGS(STRUCTURE_TO_3DI.out.mapping)
        PREP_FS_SEQS(MERGE_MAPPINGS.out.mapping)
    
    emit: 
        foldseek_db = PREP_FS_SEQS.out.foldseek_db

}