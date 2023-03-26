process  GENERATE_LIBRARY {
    container 'luisas/foldseek_tcoffee:2'
    storeDir "${params.outdir}/libraries/foldseek/${params.target_db}/$id/"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), path(structures)

    output:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree), file("${id}.library"), path(structures) , emit: library

    script:
    """
    # Identify all pairs 
    # File A B (id) 
    create_pairs.py $seqs ${id}.pairs


    # Foldseek pairwise alignment for all pairs
    while IFS="," read -r pair1 pair2; do
        foldseek easy-search \$pair1"_alphafold.pdb" \$pair2"_alphafold.pdb" \$pair1"_"\$pair2".aln" tmp -a --format-output query,qaln,target,taln

        # Convert it to fasta aligned format 
        tr "\t" "\n" < \$pair1"_"\$pair2".aln" | sed '1~2s/^/>/'  > \$pair1"_"\$pair2"_format.aln"
    done <  ${id}.pairs


    # Merge intermediate libraries
    t_coffee -in *_format.aln -lib_only -out_lib ${id}.library 
    """

}


process LIBRARY_3DI {
    container 'luisas/foldseek_tcoffee:2'
    storeDir "${params.outdir}/libraries/foldseek/${params.target_db}/$id/"
    label 'process_medium'

    input:
    tuple val(id), file(seqs), file(mapping), file(encoded_seqs)
    file(matrix)

    output:
    tuple val(id), file("${id}.library"), emit: library
    path ".command.trace", emit: metricFile

    script:
    """
    t_coffee -in $encoded_seqs -matrix $matrix -lib_only -out_lib ${id}_3di.library 
    change_library_sequences.py ${id}_3di.library $mapping ${id}.library
    """
}


process GENERATE_SEQUENCE_LIBRARY {
    container 'luisas/foldseek_tcoffee:2'
    storeDir "${params.outdir}/libraries/sequence/$id/"
    label 'process_medium'

    input:
    tuple val(id), file(seqs)

    output:
    tuple val(id), file("${id}_sequence.library"), emit: library
    path ".command.trace", emit: metricFile


    script:
    """
    t_coffee -in $seqs -lib_only -out_lib ${id}_sequence.library 
    """
}


process MERGE_LIBRARIES{
    container 'luisas/foldseek_tcoffee:2'
    storeDir "${params.outdir}/libraries/${library_id}/${params.target_db}/$id/"
    label 'process_medium'

    input:
    tuple val(id), file(library1), file(library2)
    val(library_id)

    output:
    tuple val(id), file("${id}_merged.library"), emit: library

    script:
    """
    t_coffee -in $library1 $library2 -lib_only -out_lib ${id}_merged.library 
    """
}

process GENERATE_3D_LIBRARY {
    container 'luisas/structural_regression:20'
    storeDir "${params.outdir}/libraries/$method/${params.target_db}/$id/"
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), file(seqs), file(tree), file(structures)
    val(method)

    output:
    tuple val(id), file("${id}_${method}.library"), emit: library
    path ".command.trace", emit: metricFile


    script:
    """
    # Prep templates
    for i in `awk 'sub(/^>/, "")' ${seqs}`; do
        id_pdb=`echo \$i |  sed 's./._.g'`;  echo -e ">"\$i "_P_" "\${id_pdb}" >> template_list.txt
    done

    t_coffee -in $seqs -lib_only -method $method -usetree=$tree -template_file template_list.txt -out_lib ${id}_${method}.library 
    """
}



process PRINTLIBRARY{
    container 'luisas/foldseek_tcoffee:2'
    storeDir "${params.outdir}/libraries/${library_id}/print/$id/"
    label 'process_medium'

    input:
    tuple val(id), file(msa), file(library)
    val(library_id)

    output:
    tuple val(id), file("*.html"), emit: html

    script:
    """
    t_coffee -other_pg seq_reformat -in=$msa -struc_in=$library -output=color_pdf -out=${id}.pdf 
    """
}