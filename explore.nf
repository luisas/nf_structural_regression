/*
 * enables modules
 */
nextflow.enable.dsl = 2

/*
 * defaults parameter definitions
 */


include { split_if_contains } from './modules/functions.nf'
include { MSA_IRMSD; } from './subworkflows/msa_structure_eval.nf'

refs_ch = Channel.fromPath( params.refs ).map { item -> [ item.baseName, item] }

// Collect the STRUCTURES
if (params.alphafold) {
  structures_ch = Channel.fromPath("${params.af2_db_path}/colabfold/*/*_alphafold.pdb")
                         .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0) , item.baseName.replace("_alphafold", ""), item] }
}else{
  structures_ch = Channel.fromPath("${params.experimental_structures_path}")
                         .map { item -> [split_if_contains(item.getParent().baseName, "-ref", 0) , item.baseName.replace("_header", ""), item] }
}

//refs_ch.view()
//structures_ch.view()
structures_ch = structures_ch.map{ it -> [ it[0], it[2] ]}

/*
 * main script flow
 */
workflow pipeline {

  MSA_IRMSD(refs_ch, structures_ch )

}

workflow {
  pipeline()
}

/*
 * completion handler
 */
 workflow.onComplete {

     def msg = """\
         Pipeline execution summary
         ---------------------------
         Completed at: ${workflow.complete}
         Duration    : ${workflow.duration}
         Success     : ${workflow.success}
         workDir     : ${workflow.workDir}
         exit status : ${workflow.exitStatus}
         """
         .stripIndent()

     sendMail(to: 'luisa.santus95@gmail.com', subject: 'Structural regression pipeline', body: msg)
 }
 /*
  * Extra functions
  */

 def only_first_extension(String filename){
   return filename.split("\\.").take(2).join('.')
 }
