version 1.0

task Cutadapt {
    input {
        File read1
        File? read2
        String read1output = "cut_r1.fq.gz"
        String? read2output = if defined(read2) then "cut_r2.fq.gz" else read2
        String? format
        Array[String]+? adapter
        Array[String]+? front
        Array[String]+? anywhere
        Array[String]+? adapterRead2
        Array[String]+? frontRead2
        Array[String]+? anywhereRead2
        # FIXME: default should be set at the subworkflow level, not here. Needs to wait for cromwell fix.
        Array[String]+? adapterBoth = ["AGATCGGAAGAG"]
        # contaminations = anywhereBoth
        Array[String]+? contaminations
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
        Boolean? trimPrimer
        Boolean? stripF3
        Boolean? maq
        Boolean? bwa
        Boolean? zeroCap
        Boolean? noZeroCap
        String? reportPath
        #Int compressionLevel = 1  # This only affects outputs with the .gz suffix.
        # --compression-level has a bug in 2.4 https://github.com/marcelm/cutadapt/pull/388
        #~{"--compression-level=" + compressionLevel} \
        Boolean Z = true  # equal to compressionLevel=1  # Fixme: replace once upstream is fixed.
        Int cores = 1
        Int memory = 16  # FIXME: Insane memory. Double-check if needed.
        String dockerImage = "quay.io/biocontainers/cutadapt:2.4--py37h14c3975_0"
    }

    String read2outputArg = if (defined(read2output))
        then "mkdir -p $(dirname " + read2output + ")"
        else ""

    # FIXME: This crappy overengineering can be removed once cromwell can handle subworkflow inputs correctly.
    # Some WDL magic here to set both adapters with one setting.
    # If then else's are needed to keep the variable optional and undefined
    Array[String]+? adapterForward = if (defined(adapter) || defined(adapterBoth))
                                     then select_first([adapter, adapterBoth])
                                     else adapter
    # Check if read2 is defined before applying adapters.
    Array[String]+? adapterReverse = if (defined(read2) && (defined(adapterRead2) || defined(adapterBoth)))
                                     then select_first([adapterRead2, adapterBoth])
                                     else adapterRead2

    # Same for contaminations
    Array[String]+? anywhereForward = if (defined(anywhere) || defined(contaminations))
                                      then select_first([anywhere, contaminations])
                                      else anywhere
    Array[String]+? anywhereReverse = if (defined(read2) && (defined(anywhereRead2) || defined(contaminations)))
                                      then select_first([anywhereRead2, contaminations])
                                      else anywhereRead2

    command {
        set -e
        ~{"mkdir -p $(dirname " + read1output + ")"}
        ~{read2outputArg}
        cutadapt \
        ~{"--cores=" + cores} \
        ~{true="-Z" false="" Z} \
        ~{true="-a" false="" defined(adapterForward)} ~{sep=" -a " adapterForward} \
        ~{true="-A" false="" defined(adapterReverse)} ~{sep=" -A " adapterReverse} \
        ~{true="-g" false="" defined(front)} ~{sep=" -g " front} \
        ~{true="-G" false="" defined(frontRead2)} ~{sep=" -G " frontRead2} \
        ~{true="-b" false="" defined(anywhereForward)} ~{sep=" -b " anywhereForward} \
        ~{true="-B" false="" defined(anywhereReverse)} ~{sep=" -B " anywhereReverse} \
        --output ~{read1output} ~{"--paired-output " + read2output} \
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
        ~{true="--trim-n" false="" trimN}  \
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
        File report = if defined(reportPath)
            then select_first([reportPath])
            else stdout()
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
}
