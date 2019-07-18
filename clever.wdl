version 1.0 

import "common.wdl"
import "bwa.wdl"
task Prediction {
    input {
        IndexedBamFile bamFile
        BwaIndex bwaIndex
        String outputPath        
        Int threads = 10 
        Int mem = 10
    }   
    

    command <<< 
        set -e
        mkdir -p $(dirname ~{outputPath})
        clever \
        -T ~{threads} \
        --use_mapq \
        --sorted \
        -f \
        ~{bamFile.file} \
        ~{bwaIndex.fastaFile} \
        ~{outputPath}
    >>> 

    output {
        File predictions = "~{outputPath}/predictions.vcf"
    }   
    
    runtime {
        cpu: threads
        memory: mem
        docker: "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6"
    }   

}

task Mateclever {
    input {
        File fiteredBamFile
        File indexedFiteredBamFile 
        BwaIndex bwaIndex
        File predictions
        String outputPath
        Int threads = 10 
        Int mem = 10
        Int cleverMaxDelLength = 100000
        Int maxLengthDiff= 30
        Int maxOffset = 150 
    }

    command <<<
        set -e
        mkdir -p $(dirname ~{outputPath})
        echo ~{outputPath} ~{fiteredBamFile} ~{predictions} none > predictions.list
        mateclever \
        -T ~{threads} \
        -k \
        -f \
        -M ~{cleverMaxDelLength} \
        -z ~{maxLengthDiff} \
        -o ~{maxOffset} \
        ~{bwaIndex.fastaFile} \
        predictions.list \
        ~{outputPath}
    >>>
    
    output {
        File matecleverVcf = "~{outputPath}/deletions.vcf" 
    }
    
    runtime {
        cpu: threads
        memory: mem
        docker: "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6"
    }
}
