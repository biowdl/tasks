version 1.0

# Copyright (c) 2017 Leiden University Medical Center
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

task Awk {
	input {
		File in
		String outputPrefix = "results.tsv"
		String awk
		String? sep
		Boolean compressedInput = false

		Int threads = 1
		Int timeMinutes = 10
		# Contains bwa 0.7.17 bwakit 0.7.17.dev1 and samtools 1.10.
		String dockerImage = "quay.io/biocontainers/samtools:1.21--h96c455f_1"
	}

	command {
		set -e
		mkdir -p "$(dirname ~{outputPrefix})"

		~{true="zcat" false="cat" compressedInput} ~{in} | \
			awk ~{"-F " + sep} \
				'~{awk}' \
				> ~{outputPrefix}
	}

	output {
		File out = outputPrefix
	}

	runtime {
		cpu: threads
		memory: "1GiB"
		time_minutes: timeMinutes
		partition: "short"
		slurm_partition: "short"
		docker: dockerImage
	}

	parameter_meta {
		# inputs
		in: {description: "Input (tabular) file", category: "required"}
		awk: {description: "AWK expression to transform the input", category: "required"}
		sep: {description: "Field separator used in the input file", category: "common"}
		compressedInput: {description: "Is the input compressed?", category: "common"}
		
		outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
		threads: {description: "The number of threads to use. Only used if the split input is not set.", category: "advanced"}
		memory: {description: "The amount of memory available to the job.", category: "advanced"}
		timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
		dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

		# outputs
		out: {description: "Transformed output."}
	}
}
