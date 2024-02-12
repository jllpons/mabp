#!/usr/bin/env nextflow

/*
 * Mapping and Analysis of Breakpoints in Proteins containing Segmented Domains
 *
 * Author: Joan Lluis Pons Ramon
 * Email:  <joanlluispons@gmail.com>
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
process fastasFromUniprotIDs {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path ids

    output:
    path "aa.fasta"
    path "nt.fasta"

    script:
    """
    $scripts_dir/id2aant.py $ids
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

/*
 * Get the solved residues from the PDB file.
 */
process getSolvedResidues {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path pdb_file
    val pdb

    output:
    path "${pdb}_solved_residues.json"

    script:
    """
    $scripts_dir/getSolvedRes.py ${pdb_file} ${pdb}
    """
}

/*
 * Get the HMM profile from the Pfam accession and save it in a file.
 */
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
    echo "${pdb[0..3]}" | ${scripts_dir}/pdb2UniProtID.sh
    """
}

/*
 * Get the amino acid and nucleotide sequences from the PDB code.
 * Save them in two files: `pdb_aa.fasta` and `pdb_nt.fasta`.
 * Then, add the PDB code to the header of each sequence.
 */
process pdbAaNt {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    val uniprotIdFromPdb

    output:
    path "pdb_aa.fasta"
    path "pdb_nt.fasta"

    script:
    """
    echo "${uniprotIdFromPdb}" | ${scripts_dir}/id2aant.py -p "pdb_"
    """
}

/*
 * Append the PDB amino acid and nucleotide sequences to the UniprotKB sequences.
 */
process appendPdbAaNtToFasta {

    publishDir params.outdir, mode: 'copy'

    input:
    path pdb_aa_fasta
    path aa_fasta
    path pdb_nt_fasta
    path nt_fasta
    val pdb

    output:
    path "aa_with_pdb.fasta"
    path "nt_with_pdb.fasta"

    script:
    """
    (cat $pdb_aa_fasta; echo; cat $aa_fasta) > aa_with_pdb.fasta &&
    (cat $pdb_nt_fasta; echo; cat $nt_fasta) > nt_with_pdb.fasta
    """
}

process hmmAlign {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path aa_with_pdb_fasta
    path pfam_hmm

    output:
    path "alignment.fasta"

    script:
    """
    hmmalign --outformat afa ${pfam_hmm} ${aa_with_pdb_fasta} | sed 's/\\./-/g' > alignment.fasta
    """
}

process pal2nal {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path alignment
    path nt_with_pdb_fasta

    output:
    path "codon_alignment.fasta"

    script:
    """
    pal2nal.pl ${alignment} ${nt_with_pdb_fasta} -output fasta -codontable 11 > codon_alignment.fasta
    """
}

process gard {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path codon_alignment

    output:
    path "${codon_alignment}.GARD.json"

    script:
    """
    hyphy gard ${codon_alignment} CPU=4
    """
}

/*
 * Main workflow
 */
workflow {
    ids_ch = Channel.fromPath(params.ids, checkIfExists: true)
    (aa_fasta_ch, nt_fasta_ch) = fastasFromUniprotIDs(ids_ch)

    pdb_ch = Channel.value(params.pdb)
    pdb_file_ch = getPDBfile(pdb_ch)
    solved_residues_ch = getSolvedResidues(pdb_file_ch, pdb_ch)

    uniProtIdFromPdb_ch = pdb2UniProtID(pdb_ch)
    (pdb_aa_fasta_ch, pdb_nt_fasta_ch) = pdbAaNt(uniProtIdFromPdb_ch)
    (aa_with_pdb_fasta_ch, nt_with_pdb_fasta_ch) = appendPdbAaNtToFasta(pdb_aa_fasta_ch,
                                                                        aa_fasta_ch,
                                                                        pdb_nt_fasta_ch,
                                                                        nt_fasta_ch,
                                                                        pdb_ch)

    pfam_ch = Channel.value(params.pfam)
    pfam_hmm_ch = getPfamHmm(pfam_ch)


    alignment_ch = hmmAlign(aa_with_pdb_fasta_ch, pfam_hmm_ch)
    codon_alignment_ch = pal2nal(alignment_ch, nt_with_pdb_fasta_ch)
    gard(codon_alignment_ch)
}

workflow.onComplete {
    log.info (workflow.success ? "\ncoconut oil\n" : "Oops... something went wrong :(")
}

