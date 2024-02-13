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


if (params.pdb.length() != 5) {
    log.error "PDB code must contain the chain identifier. Example: 1A2BA"
    exit 1
}

outDir = new File(params.outdir)
if (!outDir.exists()) {
    outDir.mkdirs()
}


process retrieveSequencesFromUniProt {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path ids

    output:
    path "aa.fasta"
    path "nt.fasta"

    script:
    """
    ${scripts_dir}/id2aant.py ${ids}
    """
}

process mapPDBToUniProtID {

    input:
    val pdb

    output:
    stdout

    script:
    """
    echo "${pdb[0..3]}" | ${scripts_dir}/pdb2UniProtID.sh
    """
}

process generatePDBSequencesFasta {

    conda "conda.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    val pdbMappedToUniProtID

    output:
    path "pdb_aa.fasta"
    path "pdb_nt.fasta"

    script:
    """
    echo "${pdbMappedToUniProtID}" | ${scripts_dir}/id2aant.py -p "pdb_"
    """
}

process mergePDBWithUniProtSequences {

    publishDir params.outdir, mode: 'copy'

    input:
    path pdb_aa_fasta
    path aa_fasta
    path pdb_nt_fasta
    path nt_fasta

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
 * Subworkflow 1:
 * =============
 * Responsible for retrieving and preparing all sequence data required 
 * for the analysis.
 */
workflow prepareSequenceData {

    take:
    ids
    pdb

    main:
    (aa_fasta_ch, nt_fasta_ch) = retrieveSequencesFromUniProt(ids)

    (pdb_aa_fasta_ch,
     pdb_nt_fasta_ch) = mapPDBToUniProtID(pdb) | generatePDBSequencesFasta

    (merged_aa_fasta_ch,
     merged_nt_fasta_ch) = mergePDBWithUniProtSequences(pdb_aa_fasta_ch,
                                                        aa_fasta_ch,
                                                        pdb_nt_fasta_ch,
                                                        nt_fasta_ch)

    emit:
    merged_aa_fasta_ch
    merged_nt_fasta_ch

}


process fetchPfamHMMProfile {

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

process alignSequencesToPfamHMM {

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

process createCodonAlignment {

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
 * Subworkflow 2:
 * =============
 * Focuses on aligning sequences to prepare for recombination breakpoint
 * analysis, including amino acid sequence alignment and the construction
 * of a corresponding codon alignment.
 */
workflow constructCodonAlignment {

    take:
    pfam
    merged_aa_fasta_ch
    merged_nt_fasta_ch

    main:
    pfam_hmm_ch = fetchPfamHMMProfile(pfam)
    alignment_ch = alignSequencesToPfamHMM(merged_aa_fasta_ch, pfam_hmm_ch)
    codon_alignment_ch = createCodonAlignment(alignment_ch, merged_nt_fasta_ch)

    emit:
    codon_alignment_ch

}


process detectRecombinationBreakpoints {

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
 * Subworkflow 3:
 * =============
 * Runs the GARD algorithm to detect recombination breakpoints
 * in the prepared codon alignment.
 */
workflow executeGardAnalysis {

    take:
    codon_alignment_ch
    threads_ch

    main:
    gard_analysis_ch = detectRecombinationBreakpoints(codon_alignment_ch,
                                                      threads_ch)

    emit:
    gard_analysis_ch

}


process downloadPDBStructure {

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

process extractSolvedResidues {

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


process retrieveDomainArchitectureFromCATH {

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
 * Subworkflow 4:
 * =============
 * Compares the GARD-detected recombination breakpoints against protein
 * domain architectures retrieved from CATH.
 */
workflow analyzeGardVsDomainArchitecture {

    take:
    pdb_ch
    gard_analysis_ch

    main:
    pdb_file_ch = downloadPDBStructure(pdb_ch)
    solved_residues_ch = extractSolvedResidues(pdb_file_ch, pdb_ch)

    cath_json_ch = retrieveDomainArchitectureFromCATH(pdb_ch)


    emit:
    // Tempory
    pdb_file_ch

}


/*
 * Main workflow
 * =============
 */
workflow {

    ids_ch = Channel.fromPath(params.ids, checkIfExists: true)
    pdb_ch = Channel.of(params.pdb)
    pfam_ch = Channel.of(params.pfam)
    threads_ch = Channel.of(params.threads)

    (merged_aa_fasta_ch, merged_nt_fasta_ch) = prepareSequenceData(ids_ch,
                                                                   pdb_ch)

    codon_alignment_ch = constructCodonAlignment(pfam_ch,
                                                 merged_aa_fasta_ch,
                                                 merged_nt_fasta_ch)

    gard_analysis_ch = executeGardAnalysis(codon_alignment_ch, threads_ch)

    analyzeGardVsDomainArchitecture(pdb_ch, gard_analysis_ch)
}


workflow.onComplete {
    log.info (workflow.success ? "\ncoconut oil\n" : "Oops... something went wrong :(")
}

