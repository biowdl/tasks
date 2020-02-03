version 1.0

# Copyright (c) 2018 Sequencing Analysis Support Core - Leiden University Medical Center
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

task Build {
    input {
        Boolean disableDifferenceCover = false
        File conversionTable
        File taxonomyTree
        File nameTable
        File referenceFile
        String indexBasename = "centrifuge_index"
        String outputPrefix

        Int? offrate
        Int? ftabChars
        Int? kmerCount
        File? sizeTable

        Int threads = 5
        String memory = "20G"
        String dockerImage = "quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        centrifuge-build \
        ~{"--threads " + threads} \
        ~{true="--nodc" false="" disableDifferenceCover} \
        ~{"--offrate " + offrate} \
        ~{"--ftabchars " + ftabChars} \
        ~{"--kmer-count " + kmerCount} \
        ~{"--size-table " + sizeTable} \
        ~{"--conversion-table " + conversionTable} \
        ~{"--taxonomy-tree " + taxonomyTree} \
        ~{"--name-table " + nameTable} \
        ~{referenceFile} \
        ~{outputPrefix + "/" + indexBasename}
    }

    output {
        Array[File] outputIndex = glob(outputPrefix + "/" + indexBasename + "*.cf")
    }

    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        disableDifferenceCover: {description: "Disable use of the difference-cover sample.", category: "required"}
        conversionTable: {description: "List of UIDs (unique ID) and corresponding taxonomic IDs.", category: "required"}
        taxonomyTree: {description: "Taxonomic tree (e.g. nodes.dmp).", category: "required"}
        nameTable: {description: "Name table (e.g. names.dmp).", category: "required"}
        referenceFile: {description: "A comma-separated list of FASTA files containing the reference sequences to be aligned to.", category: "required"}
        indexBasename: {description: "The basename of the index files to write.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        offrate: {description: "The number of rows marked by the indexer.", category: "common"}
        ftabChars: {description: "Calculate an initial BW range with respect to this character.", category: "common"}
        kmerCount: {description: "Use <int> as kmer-size for counting the distinct number of k-mers in the input sequences.", category: "common"}
        sizeTable: {description: "List of taxonomic IDs and lengths of the sequences belonging to the same taxonomic IDs.", category: "common"}
        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputIndex: {description: "Generated Centrifuge index."}
    }
}

task Classify {
    input {
        String inputFormat = "fastq"
        Boolean phred64 = false
        Int minHitLength = 22
        Array[File]+ indexFiles
        Array[File]+ read1
        String outputPrefix
        Array[File] read2 = []

        Int? trim5
        Int? trim3
        Int? reportMaxDistinct
        String? hostTaxIDs
        String? excludeTaxIDs

        Int threads = 4
        String memory = "16G"
        String dockerImage = "quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5"
    }

    Map[String, String] inputFormatOptions = {"fastq": "-q", "fasta": "-f", "qseq": "--qseq", "raw": "-r", "sequences": "-c"}

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        indexBasename="$(basename ~{sub(indexFiles[0], "\.[0-9]\.cf", "")})"
        for file in ~{sep=" " indexFiles}
        do
            ln ${file} $PWD/"$(basename ${file})"
        done
        centrifuge \
        ~{inputFormatOptions[inputFormat]} \
        ~{true="--phred64" false="--phred33" phred64} \
        ~{"--min-hitlen " + minHitLength} \
        ~{"--met-file " + outputPrefix + "_alignment_metrics.tsv"} \
        ~{"--threads " + threads} \
        ~{"--trim5 " + trim5} \
        ~{"--trim3 " + trim3} \
        ~{"-k " + reportMaxDistinct} \
        ~{"--host-taxids " + hostTaxIDs} \
        ~{"--exclude-taxids " + excludeTaxIDs} \
        -x $PWD/${indexBasename} \
        ~{true="-1" false="-U" length(read2) > 0} ~{sep="," read1} \
        ~{true="-2" false="" length(read2) > 0} ~{sep="," read2} \
        ~{"-S " + outputPrefix + "_classification.tsv"} \
        ~{"--report-file " + outputPrefix + "_output_report.tsv"}
    >>>

    output {
        File outputMetrics = outputPrefix + "_alignment_metrics.tsv"
        File outputClassification = outputPrefix + "_classification.tsv"
        File outputReport = outputPrefix + "_output_report.tsv"
    }

    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFormat: {description: "The format of the read file(s).", category: "required"}
        phred64: {description: "If set to true, Phred+64 encoding is used.", category: "required"}
        minHitLength: {description: "Minimum length of partial hits.", category: "required"}
        indexFiles: {description: "The files of the index for the reference genomes.", category: "required"}
        read1: {description: "List of files containing mate 1s, or unpaired reads.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        read2: {description: "List of files containing mate 2s.", category: "common"}
        trim5: {description: "Trim <int> bases from 5' (left) end of each read before alignment.", category: "common"}
        trim3: {description: "Trim <int> bases from 3' (right) end of each read before alignment.", category: "common"}
        reportMaxDistinct: {description: "It searches for at most <int> distinct, primary assignments for each read or pair.", category: "common"}
        hostTaxIDs: {description: "A comma-separated list of taxonomic IDs that will be preferred in classification procedure.", category: "common"}
        excludeTaxIDs: {description: "A comma-separated list of taxonomic IDs that will be excluded in classification procedure.", category: "common"}
        threads: {description: "The number of threads to be used.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputMetrics: {description: "File with Centrifuge metrics."}
        outputClassification: {description: "File with the classification results."}
        outputReport: {description: "File with a classification summary."}
    }
}

task Inspect {
    input {
        String printOption = "fasta"
        Array[File]+ indexFiles
        String outputPrefix

        Int? across

        String memory = "4G"
        String dockerImage = "quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5"
    }

    Map[String, String] outputOptions = {"fasta": "", "names": "--names", "summary": "--summary", "conversionTable": "--conversion-table", "taxonomyTree": "--taxonomy-tree", "nameTable": "--name-table", "sizeTable": "--size-table"}

    command <<<
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        indexBasename="$(basename ~{sub(indexFiles[0], "\.[0-9]\.cf", "")})"
        for file in ~{sep=" " indexFiles}
        do
            ln ${file} $PWD/"$(basename ${file})"
        done
        centrifuge-inspect \
        ~{outputOptions[printOption]} \
        ~{"--across " + across} \
        $PWD/${indexBasename} \
        > ~{outputPrefix + "/" + printOption}
    >>>

    output {
        File outputInspect = outputPrefix + "/" + printOption
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        printOption: {description: "The output option for inspect (fasta, summary, conversionTable, taxonomyTree, nameTable, sizeTable)", category: "required"}
        indexFiles: {description: "The files of the index for the reference genomes.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        across: {description: "When printing FASTA output, output a newline character every <int> bases.", category: "common"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputInspect: {description: "Output file according to output option."}
    }
}

task Download {
    input {
        String libraryPath
        Array[String]? domain
        String executable = "centrifuge-download"
        String? preCommand
        String? seqTaxMapPath
        String database = "refseq"
        String? assemblyLevel
        String? refseqCategory
        Array[String]? taxIds
        Boolean filterUnplaced = false
        Boolean maskLowComplexRegions = false
        Boolean downloadRnaSeqs = false
        Boolean modifyHeader = false
        Boolean downloadGiMap = false
    }

    # This will use centrifuge-download to download.
    # The bash statement at the beginning is to make sure
    # the directory for the SeqTaxMapPath exists.
    command {
        set -e -o pipefail
        ~{preCommand}
        ~{"mkdir -p $(dirname " + seqTaxMapPath + ")"}
        ~{executable} \
        -o ~{libraryPath} \
        ~{true='-d ' false='' defined(domain)}~{sep=','  domain} \
        ~{'-a "' + assemblyLevel + '"'} \
        ~{"-c " + refseqCategory} \
        ~{true='-t' false='' defined(taxIds)} '~{sep=',' taxIds}' \
        ~{true='-r' false='' downloadRnaSeqs} \
        ~{true='-u' false='' filterUnplaced} \
        ~{true='-m' false='' maskLowComplexRegions} \
        ~{true='-l' false='' modifyHeader} \
        ~{true='-g' false='' downloadGiMap} \
        ~{database} ~{">> " + seqTaxMapPath}
    }

    output {
        File seqTaxMap = "~{seqTaxMapPath}"
        File library = libraryPath
        Array[File] fastaFiles = glob(libraryPath + "/*/*.fna")
    }
 }

task DownloadTaxonomy {
    input {
        String centrifugeTaxonomyDir
        String executable = "centrifuge-download"
        String? preCommand
    }

    command {
        set -e -o pipefail
        ~{preCommand}
        ~{executable} \
        -o ~{centrifugeTaxonomyDir} \
        taxonomy
    }

    output {
        File taxonomyTree = centrifugeTaxonomyDir + "/nodes.dmp"
        File nameTable = centrifugeTaxonomyDir + "/names.dmp"
    }
 }

task Kreport {
    input {
        File centrifugeClassification
        String outputPrefix
        Array[File]+ indexFiles
        Boolean noLCA = false
        Boolean showZeros = false
        Boolean isCountTable = false

        Int? minimumScore
        Int? minimumLength

        String memory = "4G"
        String dockerImage = "quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5"
    }

    command <<< 
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        indexBasename="$(basename ~{sub(indexFiles[0], "\.[0-9]\.cf", "")})"
        for file in ~{sep=" " indexFiles}
        do
            ln ${file} $PWD/"$(basename ${file})"
        done
        centrifuge-kreport \
        -x $PWD/${indexBasename} \
        ~{true="--no-lca" false="" noLCA} \
        ~{true="--show-zeros" false="" showZeros} \
        ~{true="--is-count-table" false="" isCountTable} \
        ~{"--min-score " + minimumScore} \
        ~{"--min-length " + minimumLength} \
        ~{centrifugeClassification} \
        > ~{outputPrefix + "_kreport.tsv"}
    >>>

    output {
        File outputKreport = outputPrefix + "_kreport.tsv"
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        centrifugeClassification: {description: "File with Centrifuge classification results.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        indexFiles: {description: "The files of the index for the reference genomes.", category: "required"}
        noLCA: {description: "Do not report the LCA of multiple assignments, but report count fractions at the taxa.", category: "advanced"}
        showZeros: {description: "Show clades that have zero reads.", category: "advanced"}
        isCountTable: {description: "The format of the file is taxID<tab>COUNT.", category: "advanced"}
        minimumScore: {description: "Require a minimum score for reads to be counted.", category: "advanced"}
        minimumLength: {description: "Require a minimum alignment length to the read.", category: "advanced"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputKreport: {description: "File with kraken style report."}
    }
}

task KTimportTaxonomy {
    input {
        File inputFile
        String outputPrefix

        String memory = "4G"
        String dockerImage = "biocontainers/krona:v2.7.1_cv1"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"
        cat ~{inputFile} | cut -f 1,3 > kronaInput.krona
        ktImportTaxonomy kronaInput.krona
        cp taxonomy.krona.html ~{outputPrefix + "_krona.html"}
    }

    output {
        File outputKronaPlot = outputPrefix + "_krona.html"
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputFile: {description: "File with Centrifuge classification results.", category: "required"}
        outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputKronaPlot: {description: "Krona taxonomy plot html file."}
    }
}
