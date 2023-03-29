include { FOLDSEEK_LIBRARY; SEQUENCE_LIBRARY; TMALIGN_LIBRARY;SAP_LIBRARY } from '../subworkflows/prep_libraries.nf'
include { MERGE_LIBRARIES } from '../modules/library.nf'
include { ALIGN_WITH_LIBRARY } from '../modules/align.nf'
include { EVALUATE_MSA } from '../subworkflows/evaluate.nf'

workflow MULTI_LIB_ALIGN{

    take: 
        seqs_and_trees
        refs
        library_method
        LIBRARY_ONE
        LIBRARY_TWO
        
    main:
        // merge libraries
        libraries = LIBRARY_ONE.combine(LIBRARY_TWO, by:0)
        MERGE_LIBRARIES(libraries,library_method)

        // align
        seqs_and_trees = seqs_and_trees.map{ it -> [ it[1], it[2], it[3], it[4]]}
        seqs_trees_libraries = seqs_and_trees.combine(MERGE_LIBRARIES.out.library, by:0)

        seqs_trees_libraries.view()
         ALIGN_WITH_LIBRARY(seqs_trees_libraries, library_method)
        
        // evaluate
        if (params.evaluate){
            EVALUATE_MSA( ALIGN_WITH_LIBRARY.out.alignmentFile, refs)
        }

}