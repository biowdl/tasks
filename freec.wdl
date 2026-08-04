version 1.0

# Copyright (c) 2026 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task Freec {
    input {
        File bamFile
        File bamIndex
        File referenceFai
        File mappability
        Array[File]+ chrFiles
        String outputDir

        String sex = "XX"

        String memory = "16GiB"
        String dockerImage = "quay.io/biocontainers/control-freec:11.6b--hdbdd923_0"
        Int timeMinutes = 480
    }

    command <<<
    echo "
    [general]

    ploidy = 2,4
    sex=~{sex}

    window=100000
    step=50000

    breakPointThreshold=0.65
    breakPointType=2

    samtools=/usr/local/bin/samtools

    chrLenFile=~{referenceFai}
    chrFiles=$(dirname ~{chrFiles[0]})
    gemMappabilityFile=~{mappability}
    uniqueMatch=FALSE

    contaminationAdjustment=TRUE
    forceGCcontentNormalisation=1

    outputDir=~{outputDir}
    BedGraphOutput=FALSE

    [sample]

    mateFile=~{bamFile}
    inputFormat=BAM
    matesOrientation=FR
    " > ./freec_config

    freec -conf ./freec_conf
    >>>

    output {
        File gcProfile = "~{outputDir}/GC_profile.100000bp.cnp"
        File cnv = "~{outputDir}/~{basename(bamFile)}_CNVs"
        File info = "~{outputDir}/~{basename(bamFile)}_info.txt"
        File ratio = "~{outputDir}/~{basename(bamFile)}_ratio.txt"
        File sampleCpn= "~{outputDir}/~{basename(bamFile)}_sample.cpn"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes # !UnknownRuntimeKey
        docker: dockerImage
    }

    parameter_meta {
        bamFile: {description: "The bam file to analyse.", category: "required"}
        bamIndex: {description: "The index for the bam file.", category: "required"}
        referenceFai: {description: "The index for the reference fasta file.", category: "required"}
        mappability: {description: "The gem mappability file for the reference genome.", category: "required"}
        chrFiles: {description: "Fasta files containing one chromosome each.", category: "required"}
        outputDir: {description: "The directory to write the output to.", category: "required"}
        sex: {description: "The sample's sex.", category: "common"}
        
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}