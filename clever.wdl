version 1.0 

import "common.wdl"

task CallSV {
    input {
        IndexedBamFile bamFile
        Reference reference
        String outputPath        
        Int threads = 10 
    }   
    

    command <<< 
        clever -T ~{threads} --use_mapq --sorted -I -f ~{bamFile.file} ~{reference.fasta} ~{outputPath}
        echo ~{outputPath} ~{bamFile} ~{outputPath}/predictions.vcf} none > ~{outputPath}.list
        mateclever mateclever -k -f -M 100000 -z 30 -o 150 ~{reference} ~{outputPath}.list ~{outputPath}
    >>> 

    output {
        File cleverVcf = "~{outputPath}/predictions.vcf"
        File cleverList = "~{outputPath}.list"
        File matecleverVcf = "~{outputPath}/deletetions.vcf" 
    }   
    
    runtime {
        cpu: threads
        docker: "quay.io/biocontainers/clever-toolkit:2.4--py37hcfe0e84_5"
    }   

}
