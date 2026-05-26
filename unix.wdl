version 1.1

task fileGrep {
	input {
		File f
		File patterns
		Boolean fixedStrings = true
		String outputPrefix
		Int threads = 1
		Int timeMinutes = 10 + ceil(size(f, "GiB"))
		String dockerImage = "quay.io/biocontainers/coreutils:9.5"
	}

	command <<<
        set -e
        mkdir -p "$(dirname ~{outputPrefix})"

	grep -f ~{patterns} \
		~{true="-F" false="" fixedStrings
	} \
		~{f} > ~{outputPrefix}
    >>>

	output {
		File out = outputPrefix
	}

	runtime {
		cpu: threads
		memory: "1GiB"
		time_minutes: timeMinutes
		docker: dockerImage
	}
}

task Awk {
	input {
		File inp
		String outputPrefix
		String awk
		String? sep
		Int threads = 1
		Int timeMinutes = 10 + ceil(size(inp, "GiB"))
		# Contains bwa 0.7.17 bwakit 0.7.17.dev1 and samtools 1.10.
		String dockerImage = "quay.io/biocontainers/samtools:1.21--h96c455f_1"
	}

	command <<<
		set -e
		mkdir -p "$(dirname ~{outputPrefix})"

		cat ~{inp} | \
		    awk \
		        ~{"-F " + sep} \
		        ~{awk} \
		        > ~{outputPrefix}
	>>>

	output {
		File out = outputPrefix
	}

	runtime {
		cpu: threads
		memory: "1GiB"
		time_minutes: timeMinutes
		docker: dockerImage
	}
}
