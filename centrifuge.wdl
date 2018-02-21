# Copyright Sequencing Analysis Support Core - Leiden University Medical Center 2017
#
# Tasks from centrifuge

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
