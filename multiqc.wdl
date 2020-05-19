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

task MultiQC {
    input {
        # Use a string here so cromwell does not relocate an entire analysis directory
        Array[File] reports
        Boolean force = false
        Boolean dirs = false
        Int? dirsDepth
        Boolean fullNames = false
        String? title
        String? comment
        String? fileName
        String outDir = "."
        String? template
        String? tag
        String? ignore
        String? ignoreSamples
        File? sampleNames
        File? fileList
        Array[String]+? exclude
        Array[String]+? module
        Boolean dataDir = false
        String? dataFormat
        Boolean zipDataDir = false
        Boolean export = false
        Boolean flat = false
        Boolean interactive = true
        Boolean lint = false
        Boolean pdf = false
        Boolean megaQCUpload = false # This must be actively enabled in my opinion. The tools default is to upload.
        File? config  # A directory
        String? clConfig
    
        String memory = "4G"

        String dockerImage = "quay.io/biocontainers/multiqc:1.7--py_1"
    }

    # This is where the reports end up. It does not need to be changed by the
    # user. It is full of symbolic links, so it is not of any use to the user
    # anyway.
    String reportDir = "reports"

    # Below code requires python 3.6 or higher.
    # This makes sure all report files are in a report directory that 
    # MultiQC can investigate.
    # This creates files in report_dir / hashed_parent / file basename.
    # By hashing the parent path we make sure there are no file colissions as 
    # files from the same directory end up in the same directory, while files 
    # from other directories get their own directory. Cromwell also uses this 
    # strategy. Using python's builtin hash is unique enough for these purposes.
    
    command {
        python3 <<CODE
        import os
        from pathlib import Path 
        from typing import List

        reports: List[str] = ["~{sep='","' reports}"]
        report_dir: Path = Path("~{reportDir}")
        
        for report in reports:
            report_path = Path(report)
            hashed_parent = str(hash(str(report_path.parent)))
            new_path = report_dir / hashed_parent / report_path.name
            if not new_path.parent.exists():
                new_path.parent.mkdir(parents=True)
            os.symlink(report, str(new_path))
        CODE

        set -e
        mkdir -p ~{outDir}
        multiqc \
        ~{true="--force" false="" force} \
        ~{true="--dirs" false="" dirs} \
        ~{"--dirs-depth " + dirsDepth} \
        ~{true="--fullnames" false="" fullNames} \
        ~{"--title " + title} \
        ~{"--comment " + comment} \
        ~{"--filename " + fileName} \
        ~{"--outdir " + outDir} \
        ~{"--template " + template} \
        ~{"--tag " + tag} \
        ~{"--ignore " + ignore} \
        ~{"--ignore-samples" + ignoreSamples} \
        ~{"--sample-names " + sampleNames} \
        ~{"--file-list " + fileList} \
        ~{true="--exclude " false="" defined(exclude)}~{sep=" --exclude " exclude} \
        ~{true="--module " false="" defined(module)}~{sep=" --module " module} \
        ~{true="--data-dir" false="--no-data-dir" dataDir} \
        ~{"--data-format " + dataFormat} \
        ~{true="--zip-data-dir" false="" zipDataDir && dataDir} \
        ~{true="--export" false="" export} \
        ~{true="--flat" false="" flat} \
        ~{true="--interactive" false="" interactive} \
        ~{true="--lint" false="" lint} \
        ~{true="--pdf" false="" pdf} \
        ~{false="--no-megaqc-upload" true="" megaQCUpload} \
        ~{"--config " + config} \
        ~{"--cl-config " + clConfig } \
        ~{reportDir}
    }

    String reportFilename = if (defined(fileName))
        then sub(select_first([fileName]), "\.html$", "")
        else "multiqc"

    output {
        File multiqcReport = outDir + "/" + reportFilename + "_report.html"
        File multiqcDataDir = outDir + "/" +reportFilename + "_data"
    }

    runtime {
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        reports: {description: "Reports which multiqc should run on.", category: "required"}
        force: {description: "Equivalent to MultiQC's `--force` flag.", category: "advanced"}
        dirs: {description: "Equivalent to MultiQC's `--dirs` flag.", category: "advanced"}
        dirsDepth: {description: "Equivalent to MultiQC's `--dirs-depth` option.", category: "advanced"}
        fullNames: {description: "Equivalent to MultiQC's `--fullnames` flag.", category: "advanced"}
        title: {description: "Equivalent to MultiQC's `--title` option.", category: "advanced"}
        comment: {description: "Equivalent to MultiQC's `--comment` option.", category: "advanced"}
        fileName: {description: "Equivalent to MultiQC's `--filename` option.", category: "advanced"}
        outDir: {description: "Directory in whihc the output should be written.", category: "common"}
        template: {description: "Equivalent to MultiQC's `--template` option.", category: "advanced"}
        tag: {description: "Equivalent to MultiQC's `--tag` option.", category: "advanced"}
        ignore: {description: "Equivalent to MultiQC's `--ignore` option.", category: "advanced"}
        ignoreSamples: {description: "Equivalent to MultiQC's `--ignore-samples` option.", category: "advanced"}
        ignoreSymlinks: {description: "Equivalent to MultiQC's `--ignore-symlinks` flag.", category: "advanced"}
        sampleNames: {description: "Equivalent to MultiQC's `--sample-names` option.", category: "advanced"}
        fileList: {description: "Equivalent to MultiQC's `--file-list` option.", category: "advanced"}
        exclude: {description: "Equivalent to MultiQC's `--exclude` option.", category: "advanced"}
        module: {description: "Equivalent to MultiQC's `--module` option.", category: "advanced"}
        dataDir: {description: "Equivalent to MultiQC's `--data-dir` flag.", category: "advanced"}
        noDataDir: {description: "Equivalent to MultiQC's `--no-data-dir` flag.", category: "advanced"}
        dataFormat: {description: "Equivalent to MultiQC's `--data-format` option.", category: "advanced"}
        zipDataDir: {description: "Equivalent to MultiQC's `--zip-data-dir` flag.", category: "advanced"}
        export: {description: "Equivalent to MultiQC's `--export` flag.", category: "advanced"}
        flat: {description: "Equivalent to MultiQC's `--flat` flag.", category: "advanced"}
        interactive: {description: "Equivalent to MultiQC's `--interactive` flag.", category: "advanced"}
        lint: {description: "Equivalent to MultiQC's `--lint` flag.", category: "advanced"}
        pdf: {description: "Equivalent to MultiQC's `--pdf` flag.", category: "advanced"}
        megaQCUpload: {description: "Opposite to MultiQC's `--no-megaqc-upload` flag.", category: "advanced"}
        config: {description: "Equivalent to MultiQC's `--config` option.", category: "advanced"}
        clConfig: {description: "Equivalent to MultiQC's `--cl-config` option.", category: "advanced"}
        finished: {description: "An array of booleans that can be used to let multiqc wait on stuff.", category: "internal_use_only"}

        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.",
                      category: "advanced"}
    }

    meta {
        WDL_AID: {
            exclude: ["finished", "dependencies"]
        }
    }
}
