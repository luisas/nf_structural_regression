#!/bin/bash nextflow
params.outdir = 'results'

process GENERATE_DYNAMIC_CONFIG {
    container 'edgano/base:latest'
    tag "Config 4 Dynamic"

    input:
    val (masterAln)
    val (masterSize)
    val (slaveAln)
    val (slaveSize)

    output:
      path "${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config", emit: configFile
      tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize), emit: configValues

    script:
    """
    echo '${masterAln} ${masterSize}' > ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    echo '${slaveAln} ${slaveSize}' >> ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    """
}
