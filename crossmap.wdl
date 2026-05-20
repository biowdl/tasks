version 1.1

task Liftover {
	input {
		File vcf
		File chain
		File referenceFasta

		String outputPath = "map.vcf"
		Boolean refConsistent = false

		Int threads = 1
		String memory = "12G"
		Int timeMinutes = 60
		String dockerImage = "quay.io/biocontainers/crossmap:0.7.3--pyhdfd78af_0"
	}

	command {
		set -e
		mkdir -p "$(dirname ~{outputPath})"

		CrossMap vcf \
			~{true="--ref-consistent" false="" refConsistent} \
			~{chain} \
			~{vcf} \
			~{referenceFasta} \
			~{outputPath}
	}

	output {
		File out = outputPath
		File unmapped = outputPath + ".unmap"
	}

	runtime {
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
		container: dockerImage
	}
}

task LiftoverBed {
	input {
		File bed
		File chain

		String outputPath = "map.bed"

		Int threads = 1
		String memory = "12G"
		Int timeMinutes = 60
		String dockerImage = "quay.io/biocontainers/crossmap:0.7.3--pyhdfd78af_0"
	}

	command {
		set -e
		mkdir -p "$(dirname ~{outputPath})"

		CrossMap bed \
			~{chain} \
			~{bed} \
			~{outputPath}
	}

	output {
		File out = outputPath
		File unmapped = outputPath + ".unmap"
	}

	runtime {
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
		container: dockerImage
	}
}

