process MERGE_MAPPINGS {

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/alphabet/3di/${params.targetDB}/${id}/"
  label 'process_small'

  input:
  tuple val(id), file(files)

  output:
  tuple val(id), file("${id}.mapping"), emit: mapping

  script:
  """
  for file in $files; do cat \$file >>  "${id}.mapping"; done
  """
}


process SAVE_MERGED_DIR{
    
    storeDir "${params.outdir}/DB/${folder_name}/${db}/"

    input:
    tuple val(id), val(db), file(dirs)
    val(folder_name)

    output:
    tuple val(id), val(db), file("$id"), emit: all_files

    script:
    """
    mkdir -p $id
    for dir in $dirs; do cp \$dir/* $id/ ; done
    """
}