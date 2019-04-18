version 1.0 

import "common.wdl"

task Prediction {
    input {
        IndexedBamFile bamFile
        Reference reference
        String outputPath        
        Int threads = 10 
    }   
    

    command <<< 
        clever \
        -T ~{threads} \
        --use_mapq \
        --sorted \
        -f ~{bamFile.file} \
        ~{reference.fasta} \
        ~{outputPath}
    >>> 

    output {
        File predictions = "~{outputPath}/predictions.vcf"
    }   
    
    runtime {
        cpu: threads
        docker: "quay.io/biocontainers/clever-toolkit:2.4--py37hcfe0e84_5"
    }   

}

task Mateclever {
    input {
        IndexedBamFile bamFile
        Reference reference
        File predictions
        String outputPath
        Int threads = 10 
        Int cleverMaxDelLength = 100000
        Int maxLengthDiff= 30
        Int maxOffset = 150 
    }

    command <<<
        set -e
        echo ~{outputPath} ~{bamFile.file} ~{predictions} none > ~{outputPath}.list
        mateclever \
        -T ~{threads} \
        -k \
        -f \
        -M ~{cleverMaxDelLength} \
        -z ~{maxLengthDiff} \
        -o ~{maxOffset} \
        ~{reference.fasta} \
        ~{outputPath}.list ~{outputPath}
    >>>
    
    output {
        File cleverList = "~{outputPath}.list"
        File matecleverVcf = "~{outputPath}/deletions.vcf" 
    }
    
    runtime {
        cpu: threads
        docker: "quay.io/biocontainers/clever-toolkit:2.4--py37hcfe0e84_5"

    }
}
