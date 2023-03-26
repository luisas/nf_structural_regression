#!/bin/bash nextflow
//params.outdir = 'results'

include { set_templates_path } from './functions.nf'
path_templates = set_templates_path()


process MMSEQS_PREP_DB {
    container 'soedinglab/mmseqs2'
    storeDir "${params.dbdir}/dbs/${params.target_db}"
    label 'process_small'
    tag "$db_id"

    input:
    tuple val(db_id), file(db)

    output:
    path db

    script:
    """
    mmseqs createindex $db/$db_id tmp
    """
}


process MMSEQS_SEARCH {
    container 'luisas/mmseqs2test'
    storeDir "${params.outdir}/structures/search_hits/mmseqs/$id/${id}.${db_id}"
    label 'process_medium_high'
    tag "$id in $db_id"

    input:
    tuple val(id), val(tree_method), file(seqs), file(guide_tree)
    tuple val(db_id), file(db)

    output:
    tuple val(id), val(db_id), path("hits.m8"), emit: hits

    script:
    """
    mmseqs easy-search --min-seq-id ${params.min_id_mmseqs} -c ${params.min_cov_mmseqs} --cov-mode ${params.covmode_mmseqs} ${seqs} ${db}/${db_id} hits.m8 tmp --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov
    """

}


process TEMPLATE_FROM_DB_HITS {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/search_hits/mmseqs/$id/${id}.${db_id}"
    label 'process_low'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits)

    output:
    tuple val(id), val(db_id), file("filtered_hits.m8"), file("template.txt"), file("ids_to_download.txt"), emit: filtered_hits

    script:
    """
    python3 "${path_templates}/scripts/filter_hits.py" ${hits} "filtered_hits.m8"
    """
}


// Currently bad, no caching but just or trying out things TODO --> fix
process FETCH_STRUCTURES {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/fetched/${db_id}/${id}/"
    label 'process_small'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download)

    output:
    tuple val(id), val(db_id), file(hits), file(template), file("*_ref.pdb"), emit: fetched_structures

    script:
    """
    #bash "${path_templates}/scripts/fetch_${db_id}.sh"
    sort $ids_to_download | uniq > ids_to_download_uniq.txt

    for id in \$(cat ids_to_download_uniq.txt); do wget https://alphafold.ebi.ac.uk/files/\$id-model_v4.pdb; done

    # Horrible coding - please change if keeping
    # Generate links from name of the protein to the one matched_ref.pdb
    python3 "${path_templates}/scripts/getlinks.py" ${template} "make_links_tmp.sh" "-F1-model_v4.pdb"
    [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
    [ -f ./make_links.sh ] && bash ./make_links.sh
    [ -f ./make_links.sh ] && cat ./make_links.sh
    """
}


process FETCH_STRUCTURES_UNIPROT {
    container 'luisas/python:bio3'
    storeDir "${params.outdir}/structures/fetched/${db_id}/${id}/"
    label 'process_small'
    tag "$id in $db_id"

    input:
    tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download)

    output:
    tuple val(id), val(db_id), file(hits), file(template), file("*_ref.pdb"), emit: fetched_structures

    script:
    """
    #bash "${path_templates}/scripts/fetch_${db_id}.sh"

    function validate_url(){
      if [[ `wget -S --spider \$1  2>&1 | grep 'HTTP/1.1 200 OK'` ]]; then echo "true"; else echo "false";  fi
    }

    sort $ids_to_download | uniq > ids_to_download_uniq.txt

    for id in \$(cat ids_to_download_uniq.txt); do url="https://alphafold.ebi.ac.uk/files/AF-\$id-F1-model_v4.pdb"; if `validate_url \$url == "true"`; then wget \$url; else echo "does not exist"; fi ; done


    # Horrible coding - please change if keeping
    # Generate links from name of the protein to the one matched_ref.pdb
    python3 "${path_templates}/scripts/getlinks_uniprot.py" ${template} "make_links_tmp.sh" "pdb"
    [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
    [ -f ./make_links.sh ] && bash ./make_links.sh
    [ -f ./make_links.sh ] && cat ./make_links.sh
    """
}



process FETCH_FASTA{

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/structures/fetched/${db_id}/${id}/"
  label 'process_small'
  tag "$id in $db_id"

  input:
  tuple val(id), val(db_id), file(hits), file(template), file(ids_to_download)

  output:
  tuple val(id), val(db_id), file(hits), file(template), file("*_ref.fa"), emit: fastas

  script:
  """
  sort $ids_to_download | uniq > ids_to_download_uniq.txt

  for id in \$(cat ids_to_download_uniq.txt); do wget https://www.rcsb.org/fasta/entry/\$id -O \$id".fa"; done

  python3 "${path_templates}/scripts/getlinks.py" ${template} "make_links_tmp.sh" "fa"
  [ -f ./make_links_tmp.sh ] && tr ', ' ' ' < make_links_tmp.sh > make_links.sh
  [ -f ./make_links.sh ] && bash ./make_links.sh
  """

}


process PREP_STRUCTURES {
    container 'luisas/tcoffee_python'
    storeDir "${params.outdir}/structures/fetched_preprocessed/${id}/"
    label 'process_small'
    tag "$seq_id in $id"

    input:
    tuple val(id), val(seq_id), file(hits), file(template), file(pdb), file(fasta_ref)

    output:
    tuple val(id), val(seq_id), file(hits), file(template), file("${seq_id}_ref_ready.pdb"), emit: structures
    tuple val(id), file("${seq_id}.fa"), emit: fastas

    script:
    """
    # Extract the specific protein hit from the hits file, where all are stored
    tr "/" "_" < $hits > hitsfile
    awk '/^${seq_id}/' hitsfile > hits.txt

    # Extract chain from hits
    CHAIN=\$(awk '{ print \$2 }' hits.txt |  awk -F '[_]' '{ print \$2 }')

    # Extract PDB and make SEQRES and ATOM match
    t_coffee -other_pg extract_from_pdb -chain \$CHAIN -force -infile $pdb > ${seq_id}_fixed.pdb

    # Then pdb-seq to extract the fasta
    python3 ${path_templates}/scripts/pdb-seq.py ${seq_id}_fixed.pdb > ${seq_id}.fa

    # Extract the right chain from the fasta
    python3 ${path_templates}/scripts/extract_fasta_chain.py ${fasta_ref} \$CHAIN ${fasta_ref.baseName}_chain.fa

    # Then find-match only finds target sequence if 100% match
    # We want it like this also because its 100% btw the real fasta and the PDB
    SHIFT=\$(python3 ${path_templates}/scripts/find-match.py ${fasta_ref.baseName}_chain.fa ${seq_id}.fa)
    echo \$SHIFT

    # TODO: I need to find a solution for when some mismatch is happening - report

    # Extract the chunk from the template
    awk -v shift=\$SHIFT -v chain=\$CHAIN '{print "python3 ${path_templates}/scripts/extract_structure.py", \$9, \$10, \$1"_fixed.pdb", \$1"_ref_ready " chain, shift}' hits.txt > extract_structures.sh
    bash ./extract_structures.sh
    """
}



process GET_GDT {

  container 'luisas/tmscore'
  storeDir "${params.outdir}/structures/eval_gdt/${fam_id}/"
  label 'process_small'
  tag "$fam_id - $prot_id"

  input:
  tuple val(fam_id), val(prot_id), file(reference_pdb), file(pred_pdb)

  output:
  tuple val(fam_id), val(prot_id), file("${fam_id}_${prot_id}_GDT_report.txt"), emit: gdt_report

  script:
  """
  TMscore $reference_pdb $pred_pdb | grep "RMSD of" | awk '{print \"${fam_id},${prot_id},RMSD,\"\$6}' > ${fam_id}_${prot_id}_GDT_report.txt
  TMscore $reference_pdb $pred_pdb | grep "^TM-score" | awk '{print  \"${fam_id},${prot_id},TM-score,\"\$3}' >> ${fam_id}_${prot_id}_GDT_report.txt
  TMscore $reference_pdb $pred_pdb| grep "GDT-TS" | awk '{print  \"${fam_id},${prot_id},GDT-TS,\"\$2}' >> ${fam_id}_${prot_id}_GDT_report.txt
  """
}

process ADD_HEADER{
  container 'luisas/structural_regression:17'
  tag "${fam_name}"
  storeDir "${params.outdir}/structures/fetched_preprocessed/${fam_name}/"

  label "process_low"

  input:
  tuple val (fam_name), val(seq_id), path (af2_pdb)

  output:
	tuple val(fam_name), val(seq_id), path ("${seq_id}_header.pdb"), emit: pdb


  script:
  """
  # Add the headers
  mkdir pdbs
  t_coffee -other_pg extract_from_pdb -force -infile $af2_pdb > ${seq_id}_header.pdb
  """
}


process RENAME_PDB{
  container 'luisas/structural_regression:17'
  tag "${fam_name}"
  storeDir "${params.outdir}/structures/ready/${fam_name}/"

  label "process_low"

  input:
  tuple val (fam_name), val(seq_id), path (pdb_header)

  output:
	tuple val(fam_name), val(seq_id), path ("${seq_id}.pdb"), emit: pdbs


  script:
  """
  # Add the headers
  cp $pdb_header ${seq_id}.pdb
  """
}
