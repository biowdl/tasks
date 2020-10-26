version 1.0

# MIT License
#
# Copyright (c) 2020 Leiden University Medical Center
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

task SnpEff {
    input {
        File vcf
        File vcfIndex
        String genomeVersion
        File datadirZip
        String outputPath = "./snpeff.vcf"
        Boolean hgvs = true
        Boolean lof = true
        Boolean noDownstream = false
        Boolean noIntergenic = false
        Boolean noShiftHgvs = false
        Int? upDownStreamLen

        String memory = "50G"
        String javaXmx = "49G"
        Int timeMinutes = 60 #FIXME
        String dockerImage = "quay.io/biocontainers/snpeff:5.0--0"
    }

    command {
        set -e
        mkdir -p "$(dirname ~{outputPath})"
        unzip ~{datadirZip}
        snpEff -Xmx~{javaXmx} -XX:ParallelGCThreads=1 \
        -v \
        ~{genomeVersion} \
        -noDownload \
        -dataDir $PWD/data \
        ~{vcf} \
        ~{true="-hgvs" false="-noHgvs" hgvs} \
        ~{true="-lof" false="-noLof" lof} \
        ~{true="-no-downstream" false="" noDownstream} \
        ~{true="-no-intergenic" false="" noIntergenic} \
        ~{true="-noShiftHgvs" false="" noShiftHgvs} \
        ~{"-upDownStreamLen " + upDownStreamLen} \
        > ~{outputPath}
        rm -r $PWD/data
    }

    output {
        File outputVcf = outputPath
    }

    runtime {
        docker: dockerImage
        time_minutes: timeMinutes # !UnknownRuntimeKey
        memory: memory
    }

    parameter_meta {
        vcf: {description: "A VCF file to analyse.", category: "required"}
        vcfIndex: {description: "The index for the VCF file.", category: "required"}
        genomeVersion: {description: "The version of the genome to be used. The database for this genome must be present in the datadirZip.", category: "required"}
        datadirZip: {description: "A zip file containing the directory of databases. This zip file must contain a directory called `data`, with the database mentioned in the genomeVersion input as subdirectory.",
                     category: "required"}
        outputPath: {description: "The path to write the output to.", category: "common"}
        hgvs: {description: "Equivalent to `-hgvs` if true or `-noHgvs` if false.", category: "advanced"}
        lof: {description: "Equivalent to `-lof` if true or `-noLof` if false.", category: "advanced"}
        noDownstream: {description: "Equivalent to the `-no-downstream` flag.", category: "advanced"}
        noIntergenic: {description: "Equivalent to the `-no-intergenic` flag.", category: "advanced"}
        noShiftHgvs: {description: "Equivalent to the `-noShiftHgvs` flag.", category: "advanced"}
        upDownStreamLen: {descriptoin: "Equivalent to the `-upDownStreamLen` option.", category: "advanced"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        javaXmx: {description: "The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.",
                  category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }
}
