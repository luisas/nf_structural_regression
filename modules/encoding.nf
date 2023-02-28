process STRUCTURE_TO_3DI{
    container 'luisas/foldseek_tcoffee:2'
    tag "$id"
    storeDir "${params.outdir}/alphabet/3di/${id}/${id}.3di.fa"
    label 'process_small'
    
    input:
    tuple val(id), path (structures)
    
    output:
    tuple val (id), path ("*3di.fa"), emit: fasta_3di
    
    script:
    """
    # Convert structures to 3di
    
    for structure in *.pdb; do
        st_id=\$(echo \$structure | cut -d'.' -f1)
        foldseek structureto3didescriptor $structure \${st_id}_3di.fa
    done
    
    """
}