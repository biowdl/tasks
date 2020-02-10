version 1.0 
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
        File predictions = outputPath + "/predictions.vcf"
    }   
    
    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
    }   

    parameter_meta {
        bamFile: {description: "The bam file to process.", category: "required"}
        bamIndex: {description: "The index bam file.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
        bwaIndex: {description: "The BWA index files.", category: "required"}
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
        String memory = "15G"
        Int cleverMaxDelLength = 100000
        Int maxLengthDiff= 30
        Int maxOffset = 150 
        String dockerImage = "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
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
        File matecleverVcf = outputPath + "/deletions.vcf" 
    }
    
    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage 
    }

    parameter_meta {
        fiteredBam: {description: "The bam file where sequences less than 30bp were removed.", category: "required"}
        indexedFiteredBam: {description: "The index of the filtered bam file.", category: "required"}
        outputPath: {description: "The location the output VCF file should be written.", category: "common"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
        bwaIndex: {description: "The BWA index files.", category: "required"}
        predictions: {description: "The predicted deletions (VCF) from clever.", category: "required"}
    }
}
