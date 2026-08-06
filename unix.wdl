version 1.1

task fileGrep {
	input {
		File f
		File patterns
		Boolean fixedStrings = true
		Boolean extended = false
		Boolean onlyMatching = false
		Boolean invert = false
		Boolean exitOnError = true
		String outputPath
		Int threads = 1
		Int timeMinutes = 10 + ceil(size(f, "GiB"))
		String dockerImage = "quay.io/biocontainers/coreutils:9.5"
	}

	command <<<
	~{true="set -e" false="" exitOnError}
        mkdir -p "$(dirname ~{outputPath})"

	grep -f ~{patterns} \
		~{true="-F" false="" fixedStrings} \
		~{true="-E" false="" extended} \
		~{true="-o" false="" onlyMatching} \
		~{true="-v" false="" invert} \
		~{f} > ~{outputPath}

	~{true="" false="exit 0" exitOnError}
    >>>

	output {
		File out = outputPath
	}

	runtime {
		cpu: threads
		memory: "1GiB"
		time_minutes: timeMinutes
		docker: dockerImage
	}
}

task Grep {
	input {
		File f
		String pattern

		Boolean fixedStrings = false
		Boolean exitOnError = false
		Boolean extended = false
		Boolean onlyMatching = false
		Boolean invert = false

		String outputPath = "out.txt"
		Int threads = 1
		Int timeMinutes = 10 + ceil(size(f, "GiB"))
		String dockerImage = "quay.io/biocontainers/coreutils:9.5"
	}

	command <<<
	~{true="set -e" false="" exitOnError}
        mkdir -p "$(dirname ~{outputPath})"

	grep ~{pattern} \
		~{true="-F" false="" fixedStrings} \
		~{true="-E" false="" extended} \
		~{true="-o" false="" onlyMatching} \
		~{true="-v" false="" invert} \
		~{f} > ~{outputPath}

	~{true="" false="exit 0" exitOnError}
    >>>

	output {
		File out = outputPath
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
