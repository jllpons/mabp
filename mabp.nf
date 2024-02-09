#!/usr/bin/env nextflow

/*
 * Mapping and Analysis of Breakpoints in Proteins containing Segmented Domains
 *
 * Author: Joan Lluis Pons Ramon
 */

params.ids = "$baseDir/ids.txt"
params.pdb = ""
params.pfam = ""
params.outdir = "$baseDir/results"

scripts_dir = "$baseDir/scripts"

log.info """
M a P B   P I P E L I N E
=========================
UniProtKB IDs:    ${params.ids}
PDB code:         ${params.pdb}
Pfam accession:   ${params.pfam}
Output directory: ${params.outdir}
"""

/*
 * Create the output directory before running the pipeline.
 */
outDir = new File(params.outdir)
if (!outDir.exists()) {
    outDir.mkdirs()
}

/*
 * Check that the PDB code contains the chain identifier.
 */
if (params.pdb.length() != 5) {
    log.error "PDB code must contain the chain identifier. Example: 1A2BA"
    exit 1
}


/*
 * Get the amino acid and nucleotide sequences from the UniprotKB IDs.
 * Save them in two files: `aa.fasta` and `nt.fasta`.
 */
process fastas_from_UniprotIDs {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path ids

    output:
    path "aa.fasta"
    path "nt.fasta"

    script:
    """
    $scripts_dir/upId2AaNt.py $ids
    """
}

/*
 * Get the PDB file from the PDB code.
 */
process getPDBfile {

    publishDir params.outdir, mode: 'copy'

    input:
    val pdb

    output:
    path "${pdb[0..3]}.pdb"

    script:
    """
    wget https://files.rcsb.org/view/${pdb[0..3]}.pdb -O ${pdb[0..3]}.pdb
    """
}

process getPfamHmm {

    publishDir params.outdir, mode: 'copy'

    input:
    val pfam

    output:
    path "${pfam}.hmm"

    script:
    """
    wget https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/${pfam}?annotation=hmm -O ${pfam}.hmm.gz \
    && gunzip ${pfam}.hmm.gz
    """
}

/*
 * Get the UniprotKB ID from the PDB code.
 */
process pdb2UniProtID {

    input:
    val pdb

    output:
    stdout

    script:
    """
    echo "${pdb[0..3]}" | $scripts_dir/pdb2UniProtID.sh
    """
}

process pdbAaNt {

    publishDir params.outdir, mode: 'copy'


    input:
    val uniprotIdFromPdb

    output:
    path "pdb_aa.fasta"
    path "pdb_nt.fasta"

    script:
    """
    echo "${uniprotIdFromPdb}" | $scripts_dir/upId2AaNt.py -p "pdb_"
    """
}

workflow {
    ids_ch = Channel.fromPath(params.ids, checkIfExists: true)
    (aa_fasta_ch, nt_fasta_ch) = fastas_from_UniprotIDs(ids_ch)

    pdb_ch = Channel.value(params.pdb)
    pdb_file_ch = getPDBfile(pdb_ch)

    UniProtIdFromPdb_ch = pdb2UniProtID(pdb_ch)
    (pdb_aa_fasta_ch, pdb_nt_fasta_ch) = pdbAaNt(UniProtIdFromPdb_ch)

    pfam_ch = Channel.value(params.pfam)
    pfam_hmm_ch = getPfamHmm(pfam_ch)
}

workflow.onComplete {
    log.info (workflow.success ? "\ncoconut oil\n" : "Oops... something went wrong :(")
}

