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

task CleanSpliceJns {
    input {
        File SAMfile
        File referenceGenome
        String outputPrefix
        String outputDirPath
        File spliceJunctionAnnotation

        File? variantFile

        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        clean_splice_jns \
        ~{"--f=" + SAMfile} \
        ~{"--g=" + referenceGenome} \
        ~{"--o=" + outputDirPath + outputPrefix} \
        ~{"--s=" + spliceJunctionAnnotation} \
        ~{"--v=" + variantFile}
    }

    output {
        File outputCleanedSAM = outputDirPath + outputPrefix + "_clean.sam"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        SAMfile: "Input SAM file"
        referenceGenome: "Reference genome fasta file."
        outputPrefix: "Output file prefix."
        outputDirPath: "Output directory path."
        spliceJunctionAnnotation: "Splice junction file"
        variantFile: "VCF formatted file of variants"

        outputCleanedSAM: "Cleaned sam output file."
    }
}

task GetCorrectedSjsFromLog {
    input {
        File TElogFile
        String outputPrefix
        String outputDirPath

        Int cores = 1
        Int memory = 5
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        get_corrected_SJs_from_log \
        ~{TElogFile} \
        ~{outputDirPath + outputPrefix + ".tsv"}
    }

    output {
        File outputCorrectedSjs = outputDirPath + outputPrefix + ".tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        TElogFile: "TE log from TranscriptClean."
        outputPrefix: "Output file prefix."
        outputDirPath: "Output directory path."

        outputCorrectedSjs: "Formely noncanonical splice junctions in BED format."
    }
}

task GetSjsFromGtf {
    input {
        File GTFfile
        File genomeFile
        String outputPrefix
        String outputDirPath

        Int? minIntronSize = 21

        Int cores = 1
        Int memory = 8
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        get_SJs_from_gtf \
        ~{"--f=" + GTFfile} \
        ~{"--g=" + genomeFile} \
        ~{"--o=" + outputDirPath + outputPrefix + ".tsv"} \
        ~{"--minIntronSize=" + minIntronSize}
    }

    output {
        File outputSjsFile = outputDirPath + outputPrefix + ".tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        GTFfile: "Input GTF file"
        genomeFile: "Reference genome"
        outputPrefix: "Output file prefix."
        outputDirPath: "Output directory path."
        minIntronSize: "Minimum size of intron to consider a junction."

        outputSjsFile: "Extracted splice junctions."
    }
}

task GetTranscriptCleanStats {
    input {
        File transcriptCleanSAMfile
        String outputPrefix
        String outputDirPath

        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        get_TranscriptClean_stats \
        ~{transcriptCleanSAMfile} \
        ~{outputDirPath + outputPrefix}
    }

    output {
        File outputStatsFile = stdout()
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        transcriptCleanSAMfile: "Output SAM file from TranscriptClean"
        outputPrefix: "Output file prefix."
        outputDirPath: "Output directory path."

        outputStatsFile: "Summary stats from TranscriptClean run."
    }
}

task TranscriptClean {
    input {
        File SAMfile
        File referenceGenome
        String outputPrefix
        String outputDirPath

        File? spliceJunctionAnnotation
        File? variantFile
        Int? maxLenIndel = 5
        Int? maxSjOffset = 5
        Boolean? correctMismatches = true
        Boolean? correctIndels = true
        Boolean? correctSjs
        Boolean? dryRun = false
        Boolean? primaryOnly = false

        Int cores = 1
        Int memory = 25
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        TranscriptClean \
        ~{"-s " + SAMfile} \
        ~{"-g " + referenceGenome} \
        ~{"-o " + outputDirPath + outputPrefix} \
        ~{"-j " + spliceJunctionAnnotation} \
        ~{"-v " + variantFile} \
        ~{"--maxLenIndel=" + maxLenIndel} \
        ~{"--maxSJOffset=" + maxSjOffset} \
        ~{true="-m CORRECTMISMATCHES" false="-m false" correctMismatches} \
        ~{true="-i CORRECTINDELS" false="-i false" correctIndels} \
        ~{true="--correctSJs=CORRECTSJS" false="--correctSJs=false" correctSjs} \
        ~{true="--dryRun" false="" dryRun} \
        ~{true="--primaryOnly" false="" primaryOnly}
    }

    output {
        File outputTcFasta = outputDirPath + outputPrefix + "_clean.fa"
        File outputTcLog = outputDirPath + outputPrefix + "_clean.log"
        File outputTcSAM = outputDirPath + outputPrefix + "_clean.sam"
        File outputTcTElog = outputDirPath + outputPrefix + "_clean.TE.log"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        SAMfile: "Input SAM file containing transcripts to correct."
        referenceGenome: "Reference genome fasta file."
        outputPrefix: "Output file prefix."
        outputDirPath: "Output directory path."
        spliceJunctionAnnotation: "Splice junction file"
        maxLenIndel: "Maximum size indel to correct."
        maxSjOffset: "Maximum distance from annotated splice junction to correct."
        correctMismatches: "Set this to make TranscriptClean correct mismatches."
        correctIndels: "Set this to make TranscriptClean correct indels."
        correctSjs: "Set this to make TranscriptClean correct splice junctions."
        dryRun: "TranscriptClean will read in the data but don't do any correction."
        primaryOnly: "TranscriptClean will only output primary mappings of transcripts."

        outputTcFasta: "Fasta file containing corrected reads."
        outputTcLog: "Log file of TranscriptClean run."
        outputTcSAM: "SAM file containing corrected aligned reads."
        outputTcTElog: "TE log file of TranscriptClean run."
   }
}
