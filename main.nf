#!/usr/bin/env nextflow

/*
Mapping and Analysis of Breakpoints in Proteins containing Segmented Domains

Author: Joan Lluis Pons Ramon
*/

params.ids = "$baseDir/ids.txt"
params.outdir = "$baseDir/results"

log.info """
M a P B   P I P E L I N E
=========================
UniProtKB IDs:    ${params.ids}
Output directory: ${params.outdir}
"""

process getAAandNTfastas {
    input:
    path ids

    output:
    stdout

    script:
    """
    cat ${ids}
    """
}

workflow {
    UniprotIDs_ch = Channel.fromPath(params.ids, checkIfExists: true)
    IDsAAandNT_ch = getAAandNTfastas(UniprotIDs_ch)
    IDsAAandNT_ch.view { it }
}

workflow.onComplete {
    log.info (workflow.success ? "\ncoconut oil\n" : "Oops... something went wrong :(")
}

