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
params.threads = 1

scripts_dir = "$baseDir/scripts"

log.info """
M a P B   P I P E L I N E
=========================
UniProtKB IDs:    ${params.ids}
PDB code:         ${params.pdb}
Pfam accession:   ${params.pfam}
Output directory: ${params.outdir}
Threads for GARD: ${params.threads}
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

/*
 * Fetch CATH information about the domain architecture of the PDB protein.
 */
process fetchCathInfo {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    val pdb

    output:
    path "${pdb[0..3]}_cath_info.json"

    /*
     * Note: Here I'm using a `shell` block instead of a `script` block.
     *       This is because I want to use both shell and Nextflow variables.
     *       Shell variables are defined with `$`
     *       Nextflow variables are defined with `!{...}`
     */
    shell:
    '''
    pdb_lower=$(echo "!{pdb[0..3]}" | tr '[:upper:]' '[:lower:]')
    wget -O - "https://www.ebi.ac.uk/pdbe/api/mappings/$pdb_lower" | jq '.["'$pdb_lower'"].CATH' > "!{pdb[0..3]}_cath_info.json"
    '''
}

/*
 * Align the amino acid sequences with the Pfam HMM profile.
 */
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

/*
 * Build the codon alignment from the aligned amino acid sequences and the nucleotide sequences.
 */
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

/*
 * Run the GARD algorithm to detect recombination breakpoints in the codon alignment.
 */
process gard {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path codon_alignment
    val threads

    output:
    path "${codon_alignment}.GARD.json"

    script:
    """
    hyphy gard ${codon_alignment} CPU=${threads}
    """
}

/*
 * Main workflow
 */
workflow {
    /*
     * A) Take the list of UniprotKB IDs and generate 2 fasta files:
     *    `aa.fasta`: Amino acid sequences, obtained from the UniprotKB database.
     *    `nt.fasta`: Nucleotide sequences, obtained from the ENA database.
     *
     *    Fasta headers will appear as `>UniprotKBid_ENAid`.
     */
    ids_ch = Channel.fromPath(params.ids, checkIfExists: true)
    (aa_fasta_ch, nt_fasta_ch) = fastasFromUniprotIDs(ids_ch)

    /*
     * B) Take the PDB code and:
     *    B.1) Download the PDB file.
     *    B.2) Get the solved residues from the PDB file.
     *    B.3) Fetch CATH information about the domain architecture of the PDB protein.
     *    B.4) Get the UniprotKB ID from the PDB code.
     *        B.4.a) Generate 2 fasta files with the amino acid and nucleotide
     *               sequences from the PDB code.
     *        B.4.b) Concatenate the PDB amino acid and nucleotide sequences to
     *               the UniprotKB sequences.
     */
    pdb_ch = Channel.value(params.pdb)
    // B.1
    pdb_file_ch = getPDBfile(pdb_ch)
    // B.2
    solved_residues_ch = getSolvedResidues(pdb_file_ch, pdb_ch)
    // B.3
    cath_info_ch = fetchCathInfo(pdb_ch)
    // B.4
    uniProtIdFromPdb_ch = pdb2UniProtID(pdb_ch)
    // B.4.a
    (pdb_aa_fasta_ch, pdb_nt_fasta_ch) = pdbAaNt(uniProtIdFromPdb_ch)
    // B.4.b
    (aa_with_pdb_fasta_ch, nt_with_pdb_fasta_ch) = appendPdbAaNtToFasta(pdb_aa_fasta_ch,
                                                                        aa_fasta_ch,
                                                                        pdb_nt_fasta_ch,
                                                                        nt_fasta_ch,
                                                                        pdb_ch)

    /*
     * C) Take the Pfam accession and get the HMM profile.
     */
    pfam_ch = Channel.value(params.pfam)
    pfam_hmm_ch = getPfamHmm(pfam_ch)


    /*
     * D.1) Align the amino acid sequences with the Pfam HMM profile.
     *
     * D.2) Build the codon alignment from the aligned amino acid sequences
     *      and the nucleotide sequences.
     *
     * D.3) Run the GARD algorithm to detect recombination breakpoints in the
     *      codon alignment.
     */
    // D.1
    alignment_ch = hmmAlign(aa_with_pdb_fasta_ch, pfam_hmm_ch)
    // D.2
    codon_alignment_ch = pal2nal(alignment_ch, nt_with_pdb_fasta_ch)
    // D.3
    gard(codon_alignment_ch, params.threads)
}

workflow.onComplete {
    log.info (workflow.success ? "coconut oil\n" : "Oops... something went wrong :(")
}

