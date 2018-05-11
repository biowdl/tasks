# Copyright Sequencing Analysis Support Core - Leiden University Medical Center 2018
#
# Tasks from centrifuge
task build {

    File conversionTable
    File taxonomyTree
    File inputFasta
    String centrifugeIndexBase
    String? preCommand
    String? centrifugeBuildExecutable = "centrifuge-build"
    #Boolean? c = false
    Boolean? largeIndex = false
    Boolean? noAuto = false
    Int? bMax
    Int? bMaxDivn
    Boolean? noDiffCover = false
    Boolean? noRef = false
    Boolean? justRef = false
    Int? offRate
    Int? fTabChars
    File? nameTable
    File? sizeTable
    Int? seed
    Int? threads = 1
    Int? kmerCount

    command {
        set -e -o pipefail
        ${preCommand}
        ${"mkdir -p $(dirname " + centrifugeIndexBase + ")"}
        ${centrifugeBuildExecutable} \
        ${true='--large-index' false='' largeIndex} \
        ${true='--noauto' false='' noAuto} \
        ${'--bmax ' + bMax} \
        ${'--bmaxdivn ' + bMaxDivn} \
        ${true='--nodc' false='' noDiffCover} \
        ${true='--noref' false='' noRef} \
        ${true='--justref' false='' justRef} \
        ${'--offrate ' + offRate} \
        ${'--ftabchars ' + fTabChars} \
        ${'--name-table ' + nameTable } \
        ${'--size-table ' + sizeTable} \
        ${'--seed ' + seed} \
        ${'--kmer-count' + kmerCount} \
        ${'--threads ' + threads} \
        --conversion-table ${conversionTable} \
        --taxonomy-tree ${taxonomyTree} \
        ${inputFasta} \
        ${centrifugeIndexBase}
    }
    runtime {
        cpu: select_first([threads])
    }
}

task classify {
    String outputDir
    Boolean compressOutput = true
    String? preCommand
    String indexPrefix
    File? unpairedReads
    File read1
    File? read2
    Boolean? fastaInput
    String? outputFilePath = outputDir + "/centrifuge.out"
    String? reportFilePath = outputDir + "/centrifuge_report.tsv"
    String? metFilePath # If this is specified, the report file is empty
    Int? assignments
    Int? minHitLen
    Int? minTotalLen
    Array[String]? hostTaxIds
    Array[String]? excludeTaxIds
    Int? threads
    Int? memory

    command {
    set -e -o pipefail
    mkdir -p ${outputDir}
    ${preCommand}
    centrifuge \
    ${"-p " + threads} \
    ${"-x " + indexPrefix} \
    ${true="-f" false="" fastaInput} \
    ${true="-k " false="" defined(assignments)} ${assignments} \
    ${true="-1 " false="-U " defined(read2)} ${read1} \
    ${"-2 " + read2} \
    ${"-U " + unpairedReads} \
    ${"--report-file " + reportFilePath} \
    ${"--min-hitlen " + minHitLen} \
    ${"--min-totallen " + minTotalLen} \
    ${true="--host-taxids " false="" defined(hostTaxIds)} ${sep=',' hostTaxIds} \
    ${true="--exclude-taxids " false="" defined(excludeTaxIds)} ${sep=',' excludeTaxIds} \
    ${true="| gzip -c > " false="-S " compressOutput}${outputFilePath}${true=".gz" false="" compressOutput}
    }

    output {
        File classifiedReads = if (compressOutput) then select_first([outputFilePath + ".gz"]) else select_first([outputFilePath])
        File reportFile = select_first([reportFilePath])
    }

    runtime {
        cpu: select_first([threads, 1])
        memory: select_first([memory, 4])
    }
}

task download {
    String libraryPath
    Array[String]? domain
    String? executable = "centrifuge-download"
    String? preCommand
    String? seqTaxMapPath
    String? database = "refseq"
    String? assemblyLevel
    String? refseqCategory
    Array[String]? taxIds
    Boolean? filterUnplaced = false
    Boolean? maskLowComplexRegions = false
    Boolean? downloadRnaSeqs = false
    Boolean? modifyHeader = false
    Boolean? downloadGiMap = false

    # This will use centrifuge-download to download.
    # The bash statement at the beginning is to make sure
    # the directory for the SeqTaxMapPath exists.
    command {
        set -e -o pipefail
        ${preCommand}
        ${"mkdir -p $(dirname " + seqTaxMapPath + ")"}
        ${executable} \
        -o ${libraryPath} \
        ${true='-d ' false='' defined(domain)}${sep=','  domain} \
        ${'-a "' + assemblyLevel + '"'} \
        ${"-c " + refseqCategory} \
        ${true='-t' false='' defined(taxIds)} '${sep=',' taxIds}' \
        ${true='-r' false='' downloadRnaSeqs} \
        ${true='-u' false='' filterUnplaced} \
        ${true='-m' false='' maskLowComplexRegions} \
        ${true='-l' false='' modifyHeader} \
        ${true='-g' false='' downloadGiMap} \
        ${database} ${">> " + seqTaxMapPath}
    }
    output {
        File seqTaxMap = "${seqTaxMapPath}"
        File library = libraryPath
        Array[File] fastaFiles = glob(libraryPath + "/*/*.fna")
    }
 }

task downloadTaxonomy {
    String centrifugeTaxonomyDir
    String? executable = "centrifuge-download"
    String? preCommand
    command {
        set -e -o pipefail
        ${preCommand}
        ${executable} \
        -o ${centrifugeTaxonomyDir} \
        taxonomy
    }
    output {
        File taxonomyTree = centrifugeTaxonomyDir + "/nodes.dmp"
        File nameTable = centrifugeTaxonomyDir + "/names.dmp"
    }
 }

task kreport {
    String? preCommand
    String outputDir
    File centrifugeOut
    Boolean inputIsCompressed
    String kreportPath=sub(centrifugeOut, "\\.out.*", "\\.kreport")
    String indexPrefix
    Boolean? onlyUnique
    Boolean? showZeros
    Boolean? isCountTable
    Int? minScore
    Int? minLength
    Int? cores
    Int? memory

    command {
    set -e -o pipefail
    ${preCommand}
    centrifuge-kreport \
    -x ${indexPrefix} \
    ${true="--only-unique" false="" onlyUnique} \
    ${true="show-zeros" false="" showZeros} \
    ${true="--is-count-table" false="" isCountTable} \
    ${"--min-score " + minScore} \
    ${"--min-length " + minLength} \
    ${true="<(zcat " false="" inputIsCompressed}${centrifugeOut}${true=")" false="" inputIsCompressed} \
    > ${kreportPath}
    }

    output {
        File kreport = select_first([kreportPath])
    }

    runtime {
        cpu: select_first([cores, 1])
        memory: select_first([memory, 4])
    }
}
