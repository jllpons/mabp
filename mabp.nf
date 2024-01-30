#!/usr/bin/env nextflow

/*
 * Mapping and Analysis of Breakpoints in Proteins containing Segmented Domains
 *
 * Author: Joan Lluis Pons Ramon
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

/*
 * Get the amino acid and nucleotide sequences from the UniprotKB IDs.
 * Save them in two files: `aa.fasta` and `nt.fasta`.
 */
process fastas_from_UniprotIDs {

    conda "conda.yml"

    input:
    path ids

    output:
    path "{params.outdir}/aa.fasta"
    path "{params.outdir}/nt.fasta"

    script:
    """
    $scripts_dir/ids2aant.py $ids -o {params.outdir}
    """
}

workflow {
    ids_ch = Channel.fromPath(params.ids, checkIfExists: true)
    (aa_fasta_ch, nt_fasta_ch) = fastas_from_UniprotIDs(ids_ch)
}

workflow.onComplete {
    log.info (workflow.success ? "\ncoconut oil\n" : "Oops... something went wrong :(")
}

