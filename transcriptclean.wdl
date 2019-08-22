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

task transcriptclean {
    input {
        File SAMfile
        File referenceFile
        String outputPrefix
        String outputDirPath
        File? spliceJnsFile
        File? variantsFile
        Int? maxLenIndel
        Int? maxSJoffset
        Boolean? correctMismatches
        Boolean? correctIndels
        Boolean? correctSJs
        Boolean? dryRun
        Boolean? primaryOnly
        Int cores = 1
        Int memory = 25
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        TranscriptClean \
            ~{"-s " + SAMfile} \
            ~{"-g " + referenceFile} \
            ~{"-o " + outputDirPath + outputPrefix} \
            ~{"-j " + spliceJnsFile} \
            ~{"-v " + variantsFile} \
            ~{"--maxLenIndel=" + maxLenIndel} \
            ~{"--maxSJOffset=" + maxSJoffset} \
            ~{true="-m CORRECTMISMATCHES" false="-m false" correctMismatches} \
            ~{true="-i CORRECTINDELS" false="-i false" correctIndels} \
            ~{true="--correctSJs=CORRECTSJS" false="--correctSJs=false" correctSJs} \
            ~{true="--dryRun" false="" dryRun} \
            ~{true="--primaryOnly" false="" primaryOnly}
    }

    output {
        File outputTCfasta = outputDirPath + outputPrefix + "_clean.fa"
        File outputTClog = outputDirPath + outputPrefix + "_clean.log"
        File outputTCsam = outputDirPath + outputPrefix + "_clean.sam"
        File outputTCteLog = outputDirPath + outputPrefix + "_clean.TE.log"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }
}

task cleansplicejns {
    input{
        File SAMfile
        File referenceFile
        String outputPrefix
        String outputDirPath
        File spliceJNsFile
        File? variantsFile
        Int cores = 1
        Int memory = 4
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        clean_splice_jns \
            ~{"--f=" + SAMfile} \
            ~{"--g=" + referenceFile} \
            ~{"--o=" + outputDirPath + outputPrefix} \
            ~{"--s=" + spliceJNsFile} \
            ~{"--v=" + variantsFile}
    }

    output {
        File outputCleanSJsam = outputDirPath + outputPrefix + "_clean.sam"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }
}

task getsjsfromgtf {
    input {
        File GTFfile
        File referenceFile
        String outputPrefix
        String outputDirPath
        Int? minIntronSize
        Int cores = 1
        Int memory = 8
        String dockerImage = "biocontainers/transcriptclean:v1.0.7_cv1"
    }

    command {
        set -e pipefail
        mkdir -p ~{outputDirPath}
        get_SJs_from_gtf \
            ~{"--f=" + GTFfile} \
            ~{"--g=" + referenceFile} \
            ~{"--o=" + outputDirPath + outputPrefix + ".tsv"} \
            ~{"--minIntronSize=" + minIntronSize}
    }

    output {
        File outputSJsFile = outputDirPath + outputPrefix + ".tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }
}

task getcorrectedsjsfromlog {
    input{
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
        File outputCorrectedSJs = outputDirPath + outputPrefix + ".tsv"
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }
}

task gettranscriptcleanstats {
    input{
        File minimapSAMfile
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
            ~{minimapSAMfile} \
            ~{outputDirPath + outputPrefix}
    }

    output {
        File outputStat = stdout()
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }
}
