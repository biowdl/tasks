version 1.0

task ConvertToBin {
    input {
        File inputVcf
        String prefix = "./geva.convert"

        String memory = "4GiB"
        Int timeMinutes = 30
        Int diskGb = 1 + ceil(2.1 * size(inputVcf, "G"))
        String dockerImage = "quay.io/davycats/pkalbers-geva:5363c3db11c6b2ea2e24528affb6b68b0a939df4"
    }

    command {
        geva_v1beta \
        --vcf ~{inputVcf} \
        --out ~{prefix}
    }

    output {
        File bin = "~{prefix}.bin"
        File sample = "~{prefix}.sample.txt"
        File marker = "~{prefix}.marker.txt"
        File log = "~{prefix}.log"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        disks: "local-disk ~{diskGb} SSD" # Based on an example in dxCompiler docs
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputVcf: {description: "A VCF file (containing a single chromosome) to be converted into GEVA's binary format.", category: "required"}
        prefix: {description: "Prefix (including path) for the output files.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        diskGb: {description: "The amount of disk space needed for this job in GiB.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        bin: {description: "GEVA's binary represnetation of the VCF file."}
        sample: {description: "Sample ID mapping and metadata."}
        marker: {description: "Variant ID mapping and metadata."}
        log: {description: "GEVA's log."}
    }
}

task EstimateAge {
    input {
        File pairs
        Int effectivePopulation = 10000
        String outputPath = "./geva.estimate.txt"

        String memory = "4GiB"
        Int timeMinutes = 30
        Int diskGb = 1 + ceil(2 * size(pairs, "G"))
        String dockerImage = "quay.io/davycats/pkalbers-geva:5363c3db11c6b2ea2e24528affb6b68b0a939df4"
    }

    command {
        set -e
        mkdir ./geva_estimate_tmp
        cp ~{pairs} ./geva_estimate_tmp/tmp.pairs.txt
        Rscript /share/geva/estimate.R ./geva_estimate_tmp/tmp.pairs.txt ~{effectivePopulation}
        cp ./geva_estimate_tmp/tmp.sites2.txt ~{outputPath}
    }

    output {
        File estimates = outputPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        disks: "local-disk ~{diskGb} SSD" # Based on an example in dxCompiler docs
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        pairs: {description: "The pairs output file from GEVA.", category: "required"}
        effectivePopulation: {description: "Effective population size.", category: "advanced"}
        outputPath: {description: "Path for the output file.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        diskGb: {description: "The amount of disk space needed for this job in GiB.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        estimates: {description: "Age estimations."}
    }
}

task Geva {
    input {
        File bin
        Int position
        String prefix = "./geva.~{position}"

        String memory = "32GiB" # According to GEVA's README this is dependant on AF, so hard to figure out on the fly.
        Int timeMinutes = 120
        Int diskGb = 1 + ceil(2 * size(bin, "G"))
        String dockerImage = "quay.io/davycats/pkalbers-geva:5363c3db11c6b2ea2e24528affb6b68b0a939df4"
    }

    command {
        geva_v1beta \
        -i '~{bin}' \
        -o '~{prefix}' \
        --position '~{position}' \
        --hmm /share/geva/hmm/hmm_initial_probs.txt /share/geva/hmm/hmm_emission_probs.txt
    }

    output {
        File pairs = "~{prefix}.pairs.txt"
        File sites = "~{prefix}.sites.txt"
        File log = "~{prefix}.log"
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        disks: "local-disk ~{diskGb} SSD" # Based on an example in dxCompiler docs
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        bin: {description: "The input bin file.", category: "required"}
        position: {description: "The position to estimate the age for.", category: "required"}
        prefix: {description: "Prefix (including path) for the output files.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        diskGb: {description: "The amount of disk space needed for this job in GiB.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        pairs: {description: "Pairwise analysis results."}
        sites: {description: "Age estimations."}
        log: {description: "GEVA's log."}
    }
}
