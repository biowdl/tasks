version 1.0

# Copyright (c) 2022 Leiden University Medical Center
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

task BamReadNameToUmiTag {

    # This task processes a bam file with reads that have been extracted with
    # umi-tools extract. The UMI is extracted from the read name again and put 
    # in the bam file again with umiTag (default RX)
    input {
        File inputBam
        String outputPath = "output.bam"
        String umiTag = "RX"
        String separatorChar = "_"

        String memory = "2GiB"
        Int timeMinutes = 1 + ceil(size([inputBam], "GiB") * 10)
        String dockerImage = "quay.io/biocontainers/pysam:0.17.0--py39h051187c_0"
    }
    
    String bamIndexPath = sub(select_first([outputPath]), "\.bam$", ".bai")
    
    command <<<
        python <<CODE
        import pysam 
        import sys
        import os

        from typing import Tuple

        def split_umi_from_name(name, separator_char = "_") -> Tuple[str, str]:
            id_and_rest = name.split(maxsplit=1)
            id = id_and_rest[0]
            # If there was no whitespace id_and_rest will have length 1
            other_parts = id_and_rest[1] if len(id_and_rest) == 2 else ""
            underscore_index = id.rfind(separator_char)
            umi = id[underscore_index + 1:]
            new_id = id[:underscore_index]
            if other_parts:
                return " ".join([new_id, other_parts]), umi
            return new_id, umi

        def annotate_umis(in_file, out_file, bam_tag="RX", separator_char = "_"):
            in_bam = pysam.AlignmentFile(in_file, "rb")
            os.makedirs(os.path.dirname(out_file), exist_ok=True)
            out_bam = pysam.AlignmentFile(out_file, "wb", template=in_bam)
            # Encode bam_tag as bytes. Otherwise pysam converts it to bytes anyway.
            encoded_bam_tag = bam_tag.encode('ascii')
            for segment in in_bam:  # type: pysam.AlignedSegment
                new_name, umi = split_umi_from_name(segment.query_name, separator_char)
                segment.query_name = new_name
                # Encode umi as ascii. Otherwise pysam encodes it to bytes anyway.
                # Value type has to be a string though, otherwise pysam crashes.
                segment.set_tag(encoded_bam_tag, umi.encode('ascii'), value_type="Z")
                out_bam.write(segment)

        if __name__ == "__main__":
            annotate_umis("~{inputBam}", "~{outputPath}", "~{umiTag}", "~{separatorChar}")
            pysam.index("~{outputPath}", "~{bamIndexPath}", b=True)
        CODE
    >>>

    output {
        File outputBam = outputPath
        File outputBamIndex = bamIndexPath
    }

    runtime {
        memory: memory
        time_minutes: timeMinutes
        docker: dockerImage
    }

    parameter_meta {
        # inputs
        inputBam: {description: "The input SAM file.", category: "required"}
        outputPath: {description: "Output directory path + output file.", category: "common"}
        umiTag: {description: "The tag used for UMIs in the output BAM file.", category: "common"}
        separatorChar: {description: "Character used to separate the UMIs from the read name", category: "common"}

        memory: {description: "The amount of memory available to the job.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputBam: {description: "Sorted BAM file."}
        outputBamIndex: {description: "Sorted BAM file index."}
    }
}
