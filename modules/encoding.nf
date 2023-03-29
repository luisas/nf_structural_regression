process STRUCTURE_TO_3DI{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/alphabet/3di/${params.targetDB}/${id}/"
    label 'process_small'
    
    input:
    tuple val(id), path (structures)
    
    output:
    tuple val(id), path ("*3di.out"), emit: mapping
    path ".command.trace", emit: metricFile

    
    script:
    """
    # Convert structures to 3di
    
    for structure in *.pdb; do
        st_id=\$(echo \$structure | cut -d'.' -f1)
        foldseek structureto3didescriptor \$structure \${st_id}_3di
        cut -f1,2,3 \${st_id}_3di > \${st_id}_3di_temp.out
        echo -n `cut -f1 \${st_id}_3di_temp.out | cut -f1 -d' '` > \${st_id}_3di.out
        echo -e -n ' \t ' >> \${st_id}_3di.out
        echo -n `cut -f2 \${st_id}_3di_temp.out` >> \${st_id}_3di.out
        echo -e -n ' \t ' >> \${st_id}_3di.out
        echo -n `cut -f3 \${st_id}_3di_temp.out` >> \${st_id}_3di.out
    done
    """
}

process  ENCODE_FASTA{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/alphabet/3di_fasta/${params.targetDB}/${id}/"
    label 'process_small'
    
    input:
    tuple val(id), path(seqs), path(mapping)
    
    output:
    tuple val (id),path(seqs), path(mapping), path("${id}.3di.fa"), emit: encoded_fasta
    
    script:
    """
    encode_fasta.py $seqs $mapping ${id}.3di.fa
    """
}



process PREP_FS_SEQS{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/alphabet/fs_dir/${params.targetDB}/${id}/"
    label 'process_small'
    
    input:
    tuple val(id), path(seqs), path(mapping)
    
    output:
    tuple val (id), path("${id}_fs"), emit: fsdir
    
    script:
    """
    fs_prep.py $mapping ${id}_fs
    """
}