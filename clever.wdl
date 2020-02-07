version 1.0 

import "common.wdl"
import "bwa.wdl"
task Prediction {
    input {
        File bamFile
        File bamIndex
        BwaIndex bwaIndex
        String outputPath        
        Int threads = 10 
        String memory = "15G"
        String dockerImage = "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6"
    }   
    
    command { 
        set -e
        mkdir -p $(dirname ~{outputPath})
        clever \
        -T ~{threads} \
        --use_mapq \
        --sorted \
        -f \
        ~{bamFile} \
        ~{bwaIndex.fastaFile} \
        ~{outputPath}
    } 

    output {
        File predictions = "~{outputPath}/predictions.vcf"
    }   
    
    runtime {
        cpu: threads
        memory: mem
        docker: dockerImage
    }   

}

task Mateclever {
    input {
        File fiteredBam
        File indexedFiteredBam
        BwaIndex bwaIndex
        File predictions
        String outputPath
        Int threads = 10 
        Int mem = 15
        Int cleverMaxDelLength = 100000
        Int maxLengthDiff= 30
        Int maxOffset = 150 
        String dockerImage = "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6"
    }

    command {
        set -e
        mkdir -p $(dirname ~{outputPath})
        echo ~{outputPath} ~{fiteredBam} ~{predictions} none > predictions.list
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
    }
    
    output {
        File matecleverVcf = "~{outputPath}/deletions.vcf" 
    }
    
    runtime {
        cpu: threads
        memory: mem
        docker: dockerImage 
    }
}
