process  GENERATE_LIBRARY {
    container 'luisas/foldseek_tcoffee:2'
    storeDir "${params.outdir}/libraries/foldseek/$id/"
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