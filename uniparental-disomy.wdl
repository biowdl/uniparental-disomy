version 1.0

# Copyright (c) 2022 Sequencing Analysis Support Core - Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "tasks/gatk.wdl" as gatk
import "tasks/samtools.wdl" as samtools
import "tasks/updio.wdl" as updio

struct Trio {
    File childGvcf
    File? childGvcfIndex
    File momGvcf
    File? momGvcfIndex
    File dadGvcf
    File? dadGvcfIndex
}

workflow UniparentalDisomy {
    input {
        Array[Trio] trios
        File commonSnpBed
        File referenceFasta
        File? referenceFastaDict
        File? referenceFastaFai
        String outputDir = "."
    }

    if (!defined(referenceFastaFai) || !defined(referenceFastaDict)) {
        call samtools.DictAndFaidx as fidx {
            input:
                inputFile = referenceFasta
        }
    }
    File refFasta = select_first([fidx.outputFasta, referenceFasta])
    File refFastaDict = select_first([fidx.outputFastaDict, referenceFastaDict])
    File refFastaFai = select_first([fidx.outputFastaFai, referenceFastaFai])

    scatter (trio in trios) {
        if (!defined(trio.childGvcfIndex)) {
            call samtools.Tabix as indexChild {
                input:
                    inputFile = trio.childGvcf,
            }
        }
        if (!defined(trio.momGvcfIndex)) {
            call samtools.Tabix as indexMom {
                input:
                    inputFile = trio.momGvcf,
            }
        }
        if (!defined(trio.dadGvcfIndex)) {
            call samtools.Tabix as indexDad{
                input:
                    inputFile = trio.dadGvcf,
            }
        }
        File childGvcf = select_first([indexChild.indexedFile, trio.childGvcf])
        File momGvcf = select_first([indexMom.indexedFile, trio.momGvcf])
        File dadGvcf = select_first([indexDad.indexedFile, trio.dadGvcf])
        File childGvcfIndex = select_first([indexChild.index, trio.childGvcfIndex])
        File momGvcfIndex = select_first([indexMom.index, trio.momGvcfIndex])
        File dadGvcfIndex = select_first([indexDad.index, trio.dadGvcfIndex])

        call retrieveSamplesFromVcf as retrieveChildId {
            input:
                vcf = childGvcf
        }

        call retrieveSamplesFromVcf as retrieveMomId {
            input:
                vcf = momGvcf
        }

        call retrieveSamplesFromVcf as retrieveDadId {
            input:
                vcf = dadGvcf
        }

        String childId = retrieveChildId.samples[0]
        String momId = retrieveMomId.samples[0]
        String dadId = retrieveDadId.samples[0]
        String sampleDir = outputDir + "/" + childId

        call gatk.CombineGVCFs as CombineGVCFs {
            input:
                gvcfFiles = [childGvcf, momGvcf, dadGvcf],
                gvcfFilesIndex = [childGvcfIndex, momGvcfIndex, dadGvcfIndex],
                referenceFasta = refFasta,
                referenceFastaDict = refFastaDict,
                referenceFastaFai = refFastaFai,
                # Do not use regions here. This will make the process very slow.
        }

        call gatk.GenotypeGVCFs as GenotypeGVCFs {
            input:
                gvcfFile = CombineGVCFs.outputVcf,
                gvcfFileIndex = CombineGVCFs.outputVcfIndex,
                referenceFasta = refFasta,
                referenceFastaDict = refFastaDict,
                referenceFastaFai = refFastaFai,
                intervals = [commonSnpBed],
                # Annotation is not needed.
                annotationGroups = [],
                disableToolStandardAnnotations = true,
                outputPath = sampleDir + "/" + childId + "_trio.vcf.gz"
        }

        call updio.UpdioMultisample as runUpdio {
            input:
                multisampleVcf = GenotypeGVCFs.outputVCF,
                multisampleVcfIndex = GenotypeGVCFs.outputVCFIndex,
                childId = childId,
                momId = momId,
                dadId = dadId,
                outputPath = sampleDir
        } 
    }

    output {
        Array[File] updioFiles = flatten(runUpdio.files)
        Array[File] trioVcfs = GenotypeGVCFs.outputVCF
    }
}


task retrieveSamplesFromVcf {
    input {
        File vcf
    }

    command <<<
        bcftools query --list-samples ~{vcf}
    >>> 

    output {
        Array[String] samples = read_lines(stdout())
    }

    runtime {
        docker: "quay.io/biocontainers/bcftools:1.16--hfe4b78e_1"
        # 10 minutes should be ample as only the header and the first line are processed.
        time_minutes: 10  
        memory: "512MiB"
    }
}