version 1.0

# Copyright (c) 2017 Leiden University Medical Center
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

task Cutadapt {
    input {
        File read1
        File? read2
        String read1output = "cut_r1.fq.gz"
        String? read2output
        Array[String] adapter = []
        Array[String] front = []
        Array[String] anywhere = []
        Array[String] adapterRead2 = []
        Array[String] frontRead2 = []
        Array[String] anywhereRead2 = []
        Boolean? interleaved
        String? pairFilter
        Float? errorRate
        Boolean? noIndels
        Int? times
        Int? overlap
        Boolean? matchReadWildcards
        Boolean? noMatchAdapterWildcards
        Boolean? noTrim
        Boolean? maskAdapter
        Int? cut
        String? nextseqTrim
        String? qualityCutoff
        Int? qualityBase
        Int? length
        Boolean? trimN
        String? lengthTag
        String? stripSuffix
        String? prefix
        String? suffix
        Int? minimumLength = 2  # Necessary to prevent creation of empty reads or 1 base reads.
        Int? maximumLength
        Int? maxN
        Boolean? discardTrimmed
        Boolean? discardUntrimmed
        String? infoFilePath
        String? restFilePath
        String? wildcardFilePath
        String? tooShortOutputPath
        String? tooLongOutputPath
        String? untrimmedOutputPath
        String? tooShortPairedOutputPath
        String? tooLongPairedOutputPath
        String? untrimmedPairedOutputPath
        Boolean? colorspace
        Boolean? doubleEncode
        Boolean? stripF3
        Boolean? maq
        Boolean? bwa
        Boolean? zeroCap
        Boolean? noZeroCap
        String reportPath = "cutadapt_report.txt"
        #Int compressionLevel = 1  # This only affects outputs with the .gz suffix.
        # --compression-level has a bug in 2.4 https://github.com/marcelm/cutadapt/pull/388
        #~{"--compression-level=" + compressionLevel} \
        Boolean Z = true  # equal to compressionLevel=1  # Fixme: replace once upstream is fixed.
        Int cores = 1
        String memory = "4G"
        String dockerImage = "quay.io/biocontainers/cutadapt:2.4--py37h14c3975_0"
    }

    String realRead2output = select_first([read2output, "cut_r2.fq.gz"])
    String read2outputArg = if (defined(read2))
        then "mkdir -p $(dirname " + realRead2output + ")"
        else ""

    # FIXME: Use prefix() function for adapter, adapterRead2, etc.
    command {
        set -e
        ~{"mkdir -p $(dirname " + read1output + ")"}
        ~{read2outputArg}
        cutadapt \
        ~{"--cores=" + cores} \
        ~{true="-Z" false="" Z} \
        ~{true="-a" false="" length(adapter) > 0} ~{sep=" -a " adapter} \
        ~{true="-A" false="" length(adapterRead2) > 0} ~{sep=" -A " adapterRead2} \
        ~{true="-g" false="" length(front) > 0} ~{sep=" -g " front} \
        ~{true="-G" false="" length(frontRead2) > 0} ~{sep=" -G " frontRead2} \
        ~{true="-b" false="" length(anywhere) > 0} ~{sep=" -b " anywhere} \
        ~{true="-B" false="" length(anywhereRead2) > 0} ~{sep=" -B " anywhereRead2} \
        --output ~{read1output} ~{if defined(read2) then "-p " + realRead2output else ""} \
        ~{"--to-short-output " + tooShortOutputPath} \
        ~{"--to-short-paired-output " + tooShortPairedOutputPath} \
        ~{"--to-long-output " + tooLongOutputPath} \
        ~{"--to-long-paired-output " + tooLongPairedOutputPath} \
        ~{"--untrimmed-output " + untrimmedOutputPath} \
        ~{"--untrimmed-paired-output " + untrimmedPairedOutputPath} \
        ~{"--pair-filter " + pairFilter} \
        ~{"--error-rate " + errorRate} \
        ~{"--times " + times} \
        ~{"--overlap " + overlap} \
        ~{"--cut " + cut} \
        ~{"--nextseq-trim " + nextseqTrim} \
        ~{"--quality-cutoff " + qualityCutoff} \
        ~{"--quality-base " + qualityBase} \
        ~{"--length " + length} \
        ~{"--length-tag " + lengthTag} \
        ~{"--strip-suffix " + stripSuffix} \
        ~{"--prefix " + prefix} \
        ~{"--suffix " + suffix} \
        ~{"--minimum-length " + minimumLength} \
        ~{"--maximum-length " + maximumLength} \
        ~{"--max-n " + maxN} \
        ~{true="--discard-untrimmed" false="" discardUntrimmed} \
        ~{"--info-file " + infoFilePath } \
        ~{"--rest-file " + restFilePath } \
        ~{"--wildcard-file " + wildcardFilePath} \
        ~{true="--match-read-wildcards" false="" matchReadWildcards} \
        ~{true="--no-match-adapter-wildcards" false="" noMatchAdapterWildcards} \
        ~{true="--no-trim" false="" noTrim} \
        ~{true="--mask-adapter" false="" maskAdapter} \
        ~{true="--no-indels" false="" noIndels} \
        ~{true="--trim-n" false="" trimN} \
        ~{true="--interleaved" false="" interleaved} \
        ~{true="--discard-trimmed" false="" discardTrimmed } \
        ~{true="--colorspace" false="" colorspace} \
        ~{true="--double-encode" false="" doubleEncode} \
        ~{true="--strip-f3" false="" stripF3} \
        ~{true="--maq" false="" maq} \
        ~{true="--bwa" false="" bwa} \
        ~{true="--zero-cap" false="" zeroCap} \
        ~{true="--no-zero-cap" false="" noZeroCap} \
        ~{read1} \
        ~{read2} \
        ~{"> " + reportPath}
    }

    output{
        File cutRead1 = read1output
        File? cutRead2 = read2output
        File report = reportPath
        File? tooLongOutput=tooLongOutputPath
        File? tooShortOutput=tooShortOutputPath
        File? untrimmedOutput=untrimmedOutputPath
        File? tooLongPairedOutput=tooLongPairedOutputPath
        File? tooShortPairedOutput=tooShortPairedOutputPath
        File? untrimmedPairedOutput=untrimmedPairedOutputPath
        File? infoFile=infoFilePath
        File? restFile=restFilePath
        File? wildcardFile=wildcardFilePath
    }

    runtime {
        cpu: cores
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        read1: {
            description: "The first or single end fastq file to be run through cutadapt.",
            category: "required"
        }
        read2: {
            description: "An optional second end fastq file to be run through cutadapt.",
            category: "common"
        }
        read1output: {
            description: "The name of the resulting first or single end fastq file.",
            category: "common"
        }
        read2output: {
            description: "The name of the resulting second end fastq file.",
            category: "common"
        }
        adapter: {
            description: "A list of 3' ligated adapter sequences to be cut from the given first or single end fastq file.",
            category: "common"
        }
        front: {
            description: "A list of 5' ligated adapter sequences to be cut from the given first or single end fastq file.",
            category: "advanced"
        }
        anywhere: {
            description: "A list of 3' or 5' ligated adapter sequences to be cut from the given first or single end fastq file.",
            category: "advanced"
        }
        adapterRead2: {
            description: "A list of 3' ligated adapter sequences to be cut from the given second end fastq file.",
            category: "common"
        }
        frontRead2: {
            description: "A list of 5' ligated adapter sequences to be cut from the given second end fastq file.",
            category: "advanced"
        }
        anywhereRead2: {
            description: "A list of 3' or 5' ligated adapter sequences to be cut from the given second end fastq file.",
            category: "advanced"
        }
        interleaved: {
            description: "Equivalent to cutadapt's --interleaved flag.",
            category: "advanced"
        }
        pairFilter: {
            description: "Equivalent to cutadapt's --pair-filter option.",
            category: "advanced"
        }
        errorRate: {
            description: "Equivalent to cutadapt's --error-rate option.",
            category: "advanced"
        }
        noIndels: {
            description: "Equivalent to cutadapt's --no-indels flag.",
            category: "advanced"
        }
        times: {
            description: "Equivalent to cutadapt's --times option.",
            category: "advanced"
        }
        overlap: {
            description: "Equivalent to cutadapt's --overlap option.",
            category: "advanced"
        }
        matchReadWildcards: {
            description: "Equivalent to cutadapt's --match-read-wildcards flag.",
            category: "advanced"
        }
        noMatchAdapterWildcards: {
            description: "Equivalent to cutadapt's --no-match-adapter-wildcards flag.",
            category: "advanced"
        }
        noTrim: {
            description: "Equivalent to cutadapt's --no-trim flag.",
            category: "advanced"
        }
        maskAdapter: {
            description: "Equivalent to cutadapt's --mask-adapter flag.",
            category: "advanced"
        }
        cut: {
            description: "Equivalent to cutadapt's --cut option.",
            category: "advanced"
        }
        nextseqTrim: {
            description: "Equivalent to cutadapt's --nextseq-trim option.",
            category: "advanced"
        }
        qualityCutoff: {
            description: "Equivalent to cutadapt's --quality-cutoff option.",
            category: "advanced"
        }
        qualityBase: {
            description: "Equivalent to cutadapt's --quality-base option.",
            category: "advanced"
        }
        length: {
            description: "Equivalent to cutadapt's --length option.",
            category: "advanced"
        }
        trimN: {
            description: "Equivalent to cutadapt's --trim-n flag.",
            category: "advanced"
        }
        lengthTag: {
            description: "Equivalent to cutadapt's --length-tag option.",
            category: "advanced"
        }
        stripSuffix: {
            description: "Equivalent to cutadapt's --strip-suffix option.",
            category: "advanced"
        }
        prefix: {
            description: "Equivalent to cutadapt's --prefix option.",
            category: "advanced"
        }
        suffix: {
            description: "Equivalent to cutadapt's --suffix option.",
            category: "advanced"
        }
        minimumLength: {
            description: "Equivalent to cutadapt's --minimum-length option.",
            category: "advanced"
        }
        maximumLength: {
            description: "Equivalent to cutadapt's --maximum-length option.",
            category: "advanced"
        }
        maxN: {
            description: "Equivalent to cutadapt's --max-n option.",
            category: "advanced"
        }
        discardTrimmed: {
            description: "Equivalent to cutadapt's --quality-cutoff option.",
            category: "advanced"
        }
        discardUntrimmed: {
            description: "Equivalent to cutadapt's --discard-untrimmed option.",
            category: "advanced"
        }
        infoFilePath: {
            description: "Equivalent to cutadapt's --info-file option.",
            category: "advanced"
        }
        restFilePath: {
            description: "Equivalent to cutadapt's --rest-file option.",
            category: "advanced"
        }
        wildcardFilePath: {
            description: "Equivalent to cutadapt's --wildcard-file option.",
            category: "advanced"
        }
        tooShortOutputPath: {
            description: "Equivalent to cutadapt's --too-short-output option.",
            category: "advanced"
        }
        tooLongOutputPath: {
            description: "Equivalent to cutadapt's --too-long-output option.",
            category: "advanced"
        }
        untrimmedOutputPath: {
            description: "Equivalent to cutadapt's --untrimmed-output option.",
            category: "advanced"
        }
        tooShortPairedOutputPath: {
            description: "Equivalent to cutadapt's --too-short-paired-output option.",
            category: "advanced"
        }
        tooLongPairedOutputPath: {
            description: "Equivalent to cutadapt's --too-long-paired-output option.",
            category: "advanced"
        }
        untrimmedPairedOutputPath: {
            description: "Equivalent to cutadapt's --untrimmed-paired-output option.",
            category: "advanced"
        }
        colorspace: {
            description: "Equivalent to cutadapt's --colorspace flag.",
            category: "advanced"
        }
        doubleEncode: {
            description: "Equivalent to cutadapt's --double-encode flag.",
            category: "advanced"
        }
        stripF3: {
            description: "Equivalent to cutadapt's --strip-f3 flag.",
            category: "advanced"
        }
        maq: {
            description: "Equivalent to cutadapt's --maq flag.",
            category: "advanced"
        }
        bwa: {
            description: "Equivalent to cutadapt's --bwa flag.",
            category: "advanced"
        }
        zeroCap: {
            description: "Equivalent to cutadapt's --zero-cap flag.",
            category: "advanced"
        }
        noZeroCap: {
            description: "Equivalent to cutadapt's --no-zero-cap flag.",
            category: "advanced"
        }
        reportPath: {
            description: "The name of the file to write cutadapts's stdout to, this contains some metrics.",
            category: "common"
        }
        Z: {
            description: "Equivalent to cutadapt's -Z flag.",
            category: "advanced"
        }
        cores: {
            description: "The number of cores to use.",
            category: "advanced"
        }
        memory: {
            description: "The amount of memory this job will use.",
            category: "advanced"
        }
        dockerImage: {
            description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
            category: "advanced"
        }
    }
}
