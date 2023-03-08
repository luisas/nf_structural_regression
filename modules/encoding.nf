process STRUCTURE_TO_3DI{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/alphabet/3di/${id}/"
    label 'process_small'
    
    input:
    tuple val(id), path (structures)
    
    output:
    tuple val(id), path ("*3di.out"), emit: mapping
    
    script:
    """
    # Convert structures to 3di
    
    for structure in *.pdb; do
        st_id=\$(echo \$structure | cut -d'.' -f1)
        foldseek structureto3didescriptor \$structure \${st_id}_3di
        cut -f1,2,3 \${st_id}_3di > \${st_id}_3di.out
    done
    """
}

process  ENCODE_FASTA{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/alphabet/3di_fasta/${id}/"
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