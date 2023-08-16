#!/bin/bash nextflow

nxf_min_ver = "20.07.1"

def set_templates_path () {
    if( !nextflow.version.matches(">=${nxf_min_ver}") ) {
        println "It is advisable to run this workflow with Nextflow version $nxf_min_ver or greater -- You are running version $nextflow.version"
        path_templates = workflow.commandLine.contains('--pipeline') ? "$baseDir/modules/regressive_alignment/modules/templates" : "$baseDir/modules/templates"
    }
    else {
        path_templates = "${moduleDir}/templates"
    }
    return path_templates
}


def split_if_contains(s,sep, index ){

  if(s.contains(sep)){
    return s.split(sep)[index]
  }else{
    return s
  }

}

def static String remove_suffix( String self ) {
    if(self.contains("_")){
      return self.split("_")[0]
    }else{
      return self
    }

}
