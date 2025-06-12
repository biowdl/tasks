version 1.0

task Sniffles {
	input {
		File bam
		File bamIndex
		File? referenceFasta
		File? tandemRepeats
		File? genotypeVcfs

		String outputVcf = "output.vcf"
		String options = " "

		Boolean mosaic = false
		Boolean outputRnames = false

		Int threads = 4
		String memory = "32GiB" # 16 was fine without --no-qc
		Int timeMinutes = 1 + ceil(size(bam, "GiB") * 2)
		String dockerImage = "quay.io/biocontainers/sniffles:2.6.2--pyhdfd78af_0"
	}

	meta {
		citation: {
			doi: "10.1038/s41587-023-02024-y"
		}
	}

	command {
		set -e
		mkdir -p "$(dirname ~{outputVcf})"

		sniffles \
			~{"--reference " + referenceFasta} \
			~{"--tandem-repeats " + tandemRepeats} \
			~{"--genotype-vcf " + genotypeVcfs} \
			--threads ~{threads} \
			~{true="--output-rnames" false="" outputRnames} \
			~{true="--mosaic" false="" mosaic} \
			--input ~{bam} \
			~{options} \
			--vcf ~{outputVcf}
	}

	output {
		File vcf = outputVcf
	}

	runtime {
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
		docker: dockerImage
	}
}
