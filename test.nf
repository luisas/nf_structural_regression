#!/usr/bin/env nextflow
nextflow.enable.dsl=2

extracted_sequences = Channel.fromPath( '~/Desktop/test/sample.fa' ) \
                             .splitFasta( record: [id: true, seqString: true ]) \
                             .view()

available_structures = Channel.fromPath( '~/Desktop/test/structures/' ) \
                              .map{it -> it.baseName}.view()
