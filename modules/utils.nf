process MERGE_MAPPINGS {

  container 'luisas/python:bio3'
  storeDir "${params.outdir}/alphabet/3di/${id}/"
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