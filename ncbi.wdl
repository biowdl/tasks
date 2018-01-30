task genomeDownload {
    String outputPath
    String? section = "refseq"
    String? format = "all"
    String? assemblyLevel = "all"
    String? taxId
    String? refseqCategory
    Boolean? humanReadable
    String? ncbiBaseUri
    Int? parallel
    Int? retries
    Boolean? verbose=true
    Boolean? debug
    String? domain = "all"

    String? executable = "ncbi-genome-download"
    String? preCommand

    command {
        set -e -o pipefail
        ${preCommand}
        ${executable} \
        ${"--section " + section} \
        ${"--format " + format} \
        ${"--assembly-level " + assemblyLevel } \
        ${"--taxid " + taxId } \
        ${"--refseq-category " + refseqCategory} \
        ${"--output-folder " + outputPath } \
        ${true="--human-readable" false="" humanReadable} \
        ${"--uri " + ncbiBaseUri } \
        ${"--parallel " + parallel } \
        ${"--retries " + retries } \
        ${true="--verbose" false="" verbose } \
        ${true="--debug" false ="" debug } \
        ${domain}

        # Check md5sums for all downloaded files
        for folder in $(realpath ${outputPath})/*/*/*
            do
                (
                md5sums="$(
                    cd $folder
                    for file in *
                    do
                        if [[ ! $file == "MD5SUMS" ]]
                        then
                            grep $file MD5SUMS
                        fi
                    done
                    )"
                cd $folder; echo $md5sums | md5sum -c)
            done
    }

    output {
        Array[File] fastaGzFiles = glob(outputPath + "/*/*/*/*_genomic.fna.gz")
        Array[File] genbankGzFiles = glob(outputPath + "/*/*/*/*_genomic.gbff.gz")
        Array[File] featuresGzFiles = glob(outputPath + "/*/*/*/*_feature_table.txt.gz")
        Array[File] gffGzFiles = glob(outputPath + "/*/*/*/*_genomic.gff.gz")
        Array[File] proteinFastaGzFiles = glob(outputPath + "/*/*/*/*_protein.faa.gz")
        Array[File] genpeptGzFiles = glob(outputPath + "/*/*/*/*_protein.gpff.gz")
        Array[File] wgsGzFiles = glob(outputPath + "/*/*/*/*_wgsmaster.gbff.gz")
        Array[File] cdsFastaGzFiles = glob(outputPath + "/*/*/*/*_cds_from_genomic.fna.gz")
        Array[File] rnaFastaGzFiles = glob(outputPath + "/*/*/*/*_rna_from_genomic.fna.gz")
        Array[File] assemblyReportFiles = glob(outputPath + "/*/*/*/*_assembly_report.txt")
        Array[File] assemblyStatsFiles = glob(outputPath + "/*/*/*/*_assembly_stats.txt")
    }
 }


task downloadNtFasta{
    String libraryPath
    String seqTaxMapPath
    Boolean? unzip = true
    String ntDir = libraryPath + "/nt"
    String ntFilePath = ntDir + "/nt.fna"
    command {
        set -e -o pipefail
        mkdir -p ${ntDir}
        rsync -av --partial rsync://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz* ${ntDir}
        (cd ${ntDir} && md5sum -c nt.gz.md5)
        # Only unzip when necessary
        if ${true='true' false='false' unzip}
        then
            zcat ${ntDir}/nt.gz > ${ntFilePath}
        fi
        }
    output {
        File ntFileGz = ntDir + "/nt.gz"
        File library = libraryPath
        # Added array file to allow for multiple downloads later.
        # Also allows for easier pipeline logic.
        Array[File] ntFastas = glob(ntDir + "/*.fna")
        Array[File] ntFastasGz = glob(ntDir + "/nt*.gz")
    }
}

task downloadAccessionToTaxId {
    String downloadDir
    Boolean gzip = false
    command {
        set -e -o pipefail
        mkdir -p ${downloadDir}
        rsync -av --partial rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_*.accession2taxid.gz* ${downloadDir}
        (cd ${downloadDir} && md5sum -c *.md5)
        for file in ${downloadDir}/nucl_*.accession2taxid.gz
        do
            zcat $file | tail -n +2 | cut -f 2,3 ${true="| gzip " false='' gzip}> $file.seqtaxmap${true='.gz' false='' gzip}
        done
        }
    output {
        Array[File] seqTaxMaps = glob(downloadDir + "/*.seqtaxmap")
        Array[File] seqTaxMapsGz = glob(downloadDir + "/*.seqtaxmap.gz")
    }
}
