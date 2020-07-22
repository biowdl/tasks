version 1.0

# Copyright (c) 2019 Sequencing Analysis Support Core - Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task GetSJsFromGtf {
    input {
        File gtfFile
        File genomeFile
        String outputPrefix
        Int minIntronSize = 21

        String memory = "8G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/transcriptclean:v2.0.2_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        get_SJs_from_gtf \
        --f=~{gtfFile} \
        --g=~{genomeFile} \
        --minIntronSize=~{minIntronSize} \
        ~{"--o=" + outputPrefix + ".tsv"}
    }

    output {
        File spliceJunctionFile = outputPrefix + ".tsv"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        gtfFile: {description: "Input gtf file", category: "required"}
        genomeFile: {description: "Reference genome", category: "required"}
        minIntronSize: {description: "Minimum size of intron to consider a junction.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}
        # outputs
        spliceJunctionFile: {description: "Extracted splice junctions."}
    }
}

task GetTranscriptCleanStats {
    input {
        File inputSam
        String outputPrefix

        String memory = "4G"
        Int timeMinutes = 30
        String dockerImage = "biocontainers/transcriptclean:v2.0.2_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        get_TranscriptClean_stats \
        ~{inputSam} \
        ~{outputPrefix}
    }

    output {
        File statsFile = stdout()
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputSam: {description: "Output sam file from transcriptclean", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        statsFile: {description: "Summary stats from transcriptclean run."}
    }
}

task TranscriptClean {
    input {
        File inputSam
        File referenceGenome
        Int maxLenIndel = 5
        Int maxSJOffset = 5
        String outputPrefix
        Boolean correctMismatches = true
        Boolean correctIndels = true
        Boolean correctSJs = true
        Boolean dryRun = false
        Boolean primaryOnly = false
        Boolean canonOnly = false
        Int bufferSize = 100
        Boolean deleteTmp = true

        File? spliceJunctionAnnotation
        File? variantFile

        Int cores = 1
        String memory = "25G"
        Int timeMinutes = 2880
        String dockerImage = "biocontainers/transcriptclean:v2.0.2_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        TranscriptClean \
        -s ~{inputSam} \
        -g ~{referenceGenome} \
        -t ~{cores} \
        --maxLenIndel=~{maxLenIndel} \
        --maxSJOffset=~{maxSJoffset} \
        -o ~{outputPrefix} \
        ~{true="-m true" false="-m false" correctMismatches} \
        ~{true="-i true" false="-i false" correctIndels} \
        ~{true="--correctSJs=true" false="--correctSJs=false" correctSJs} \
        ~{true="--dryRun" false="" dryRun} \
        ~{true="--primaryOnly" false="" primaryOnly} \
        ~{true="--canonOnly" false="" canonOnly} \
        --bufferSize=~{bufferSize} \
        ~{true="--deleteTmp" false="" deleteTmp} \
        ~{"-j " + spliceJunctionAnnotation} \
        ~{"-v " + variantFile}
    }

    output {
        File fastaFile = outputPrefix + "_clean.fa"
        File logFile = outputPrefix + "_clean.log"
        File outputSam = outputPrefix + "_clean.sam"
        File logFileTE = outputPrefix + "_clean.TE.log"
    }

    runtime {
        cpu: cores
        memory: memory
        time_minute: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputSam: {description: "Input sam file containing transcripts to correct.", category: "required"}
        referenceGenome: {description: "Reference genome fasta file.", category: "required"}
        maxLenIndel: {description: "Maximum size indel to correct.", category: "advanced"}
        maxSJOffset: {description: "Maximum distance from annotated splice junction to correct.", category: "advanced"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        correctMismatches: {description: "Set this to make transcriptclean correct mismatches.", category: "common"}
        correctIndels: {description: "Set this to make transcriptclean correct indels.", category: "common"}
        correctSJs: {description: "Set this to make transcriptclean correct splice junctions.", category: "common"}
        dryRun: {description: "Transcriptclean will read in the data but don't do any correction.", category: "advanced"}
        primaryOnly: {description: "Only output primary mappings of transcripts.", category: "advanced"}
        canonOnly: {description: "Only output canonical transcripts and transcript containing annotated noncanonical junctions.", category: "advanced"}
        bufferSize: {description: "Number of lines to output to file at once by each thread during run.", category: "common"}
        deleteTmp: {description: "The temporary directory generated by transcriptclean will be removed.", category: "common"}
        spliceJunctionAnnotation: {description: "Splice junction file.", category: "common"}
        variantFile: {description: "Vcf formatted file of variants.", category: "common"}
        cores: {description: "The number of cores to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
        
        # outputs
        fastaFile: {description: "Fasta file containing corrected reads."}
        logFile: {description: "Log file of transcriptclean run."}
        outputSam: {description: "Sam file containing corrected aligned reads."}
        logFileTE: {description: "TE log file of transcriptclean run."}
   }
}
