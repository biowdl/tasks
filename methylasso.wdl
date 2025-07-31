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

task MethyLasso {
	input {
		String control_name = "control"
		Array[File] control
		String condition_name
		Array[File] condition

		String outputPrefix = "."

		Int threads = 8
		String memory = "180G"
		Int timeMinutes = 6 * 60
		String dockerImage = "ghcr.io/biowdl/docker-methylasso:sha256-e08c1bedbd0a887e8a76035bc713991e03364dcec9e8b83d7732be55f84dafeb.sig"
	}

	# Probably will work once we filter for just chr1-22
	command {
		set -e
		mkdir -p '~{outputPrefix}'

		Rscript /methylasso-main/MethyLasso.R \
			--n1 ~{control_name} \
			--c1 ~{sep="," control} \
			--n2 ~{condition_name} \
			--c2 ~{sep="," condition} \
			--threads ~{threads} \
			--meth 4 \
			--cov 5 \
			-o '~{outputPrefix}'
	}

	output {
		Array[File] comparison_table = glob("~{outputPrefix}/*vs*.tsv")
		Array[File] comparison_plot = glob("~{outputPrefix}/*vs*.pdf")
		Array[File] plots = glob("~{outputPrefix}/*.pdf")
		Array[File] tables = glob("~{outputPrefix}/*.tsv")
	}

	runtime {
		cpu: threads
		memory: memory
		time_minutes: timeMinutes
		docker: dockerImage
	}

	parameter_meta {
		# inputs
		control_name: {description: "The name of the control condition.", category: "advanced"}
		control: {description: "The set of tabular methylation levels in the control condition", category: "required"}
		condition_name: {description: "The name of the condition being explored.", category: "required"}
		condition: {description: "The set of tabular methylation levels for the experimental condition.", category: "required"}
		
		outputPrefix: {description: "Output directory path + output file prefix.", category: "required"}
		threads: {description: "The number of threads to use. Only used if the split input is not set.", category: "advanced"}
		memory: {description: "The amount of memory available to the job.", category: "advanced"}
		timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
		dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

		# outputs
		comparison_table: {description: "Set of comparison tables for results, the main tabular results you're interested in."}
		comparison_plot: {description: "PDF formatted plot of the condition vs control"}
		plots: {description: "All produced plots, showing mostly methylation level distributions"}
		tables: {description: "DMR regions: UMRs, LMRs, etc."}
	}
}
