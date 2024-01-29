#!/usr/bin/env nextflow

/*
Mapping and Analysis of Breakpoints in Proteins containing Segmented Domains

Author: Joan Lluis Pons Ramon
*/

params.ids = "$baseDir/ids.txt"
params.outdir = "$baseDir/results"

scripts_dir = "$baseDir/scripts"

log.info """
M a P B   P I P E L I N E
=========================
UniProtKB IDs:    ${params.ids}
Output directory: ${params.outdir}
"""

process AaNtFastas_from_UniprotIDs {
    // Get the amino acid and nucleotide sequences from the UniprotKB IDs.
    // Save them in two files: `aa.fasta` and `nt.fasta`.
    // Fasta headers will appear as: `>UniprotAccession_ENAAccession`
    // on both files.

    conda "conda.yml"

    input:
    path ids

    output:
    stdout

    script:
    """
    $scripts_dir/ids2aant.py -h
    """
}

workflow {
    UniprotIDs_ch = Channel.fromPath(params.ids, checkIfExists: true)
    AaNtFastas_from_UniprotIDs_ch = AaNtFastas_from_UniprotIDs(UniprotIDs_ch)
    AaNtFastas_from_UniprotIDs_ch.view { it }
}

workflow.onComplete {
    log.info (workflow.success ? "\ncoconut oil\n" : "Oops... something went wrong :(")
}

