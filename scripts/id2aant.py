#!/usr/bin/env python3

"""
From a list of UniProt IDs, this script fetches the amino acid sequences from
UniProt and the nucleotide sequences from ENA. Then, it writes them in two
fasta files. FASTA headers will appear as: `>UniprotAccession_ENAAccession`.
Accepts IDs list from sdin or from file if -f is used.

Usage:  id2aant.py <FILE> [options] or cat <FILE> | id2aant.py [options]

Author: Joan Lluis Pons Ramon
Email:  <joanlluispon@gmail.com>
"""

import argparse
import json
import os
import sys

import requests
from requests.adapters import HTTPAdapter, Retry

from API_URL import UNIPROT_JSON_ENTRY, ENA_NUCLEOTIDE_FASTA_SEQUENCE


__version__ = "0.1.1"


# <https://stackoverflow.com/questions/18275023/dont-show-long-options-twice-in-print-help-from-argparse>
class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super().__init__(prog, max_help_position=40)

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


class Protein:
    """
    Protein object with the data fetched from UniProt or ENA.
    """

    def __init__(self, input_id: str):
        self.input_id = input_id
        self.uniprot_accession = None
        self.uniprot_aa_sequence = None
        self.ena_accession = None
        self.ena_nucleotide_sequence = None

    def generate_fasta_header(self) -> str:
        """
        Generate FASTA header.

        Args:
            None

        Returns:
            header: FASTA header.
        """
        return f">{self.uniprot_accession}_{self.ena_accession}"

    def is_valid(self) -> bool:
        """
        Check if the protein object contains all the data needed to generate both
        amino acid and nucleotide FASTA sequences.

        Args:
            None

        Returns:
            valid: True if the protein object is valid, False otherwise.
        """
        if self.uniprot_accession and self.uniprot_aa_sequence and self.ena_accession and self.ena_nucleotide_sequence:
            return True
        return False


def eprint(*args, **kwargs):
    """Print to stderr"""
    print(*args, file=sys.stderr, **kwargs)


def setup_argparse() -> argparse.ArgumentParser:
    """
    Customize and setup argparse.

    Args:
        None

    Returns:
        parser: argparse.ArgumentParser object
    """

    fmt = lambda prog: CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(
            prog="id2aant",
            formatter_class=fmt,
            usage="id2aant.py <FILE> [options] or cat <FILE> | id2aant.py [options]",
            description="""
Given a list of UniProt IDs, fetches the amino acid sequences from UniProt and
the nucleotide sequences from ENA. Then, it writes them in two fasta files.
FASTA headers will appear as: `>UniprotAccession_ENAAccession`.
            """,
            epilog="""
Examples:
    id2aant.py list.txt -o results
    cat list.txt | id2aant.py
    id2aant.py list.txt -o results -p prefix

For more information, visit <https://github.com/jllpons/mabp>.
            """,
            )

    parser.add_argument(
            "file",
            metavar="FILE",
            nargs="?",
            type=str,
            help="UniProt IDs separated by newlines. Reads a file or from stdin.",
            )
    parser.add_argument(
            "-o",
            "--output",
            metavar="DIR",
            type=str,
            default=f"{os.getcwd()}",
            help="Output directory. Created if it does not exist. Default: CWD.",
            )
    parser.add_argument(
            "-p",
            "--prefix",
            metavar="STR",
            type=str,
            default="",
            help="Specify a prefix for the output files.",
            )
    parser.add_argument(
            "-q",
            "--quiet",
            action="store_true",
            default=False,
            help="Do not print anything to stderr.",
            )
    parser.add_argument(
            "-V",
            "--version",
            action="version",
            version=f"%(prog)s {__version__}",
            )

    return parser


def read_ids(file: str) -> list:
    """
    Read IDs from file.

    Args:
        file: Path to file containing IDs separated by newlines.

    Returns:
        ids: List of IDs.

    Raises:
        ValueError: If no IDs are found.
    """

    with open(file, "r") as f:
        ids = f.read().splitlines()

    if not ids:
        raise ValueError("[id2aant.py] Error: No IDs found.")

    return ids


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate arguments.

    Args:
        args: argparse.Namespace object

    Returns:
        None

    Raises:
        ValueError: If no input is provided.
        FileNotFoundError: If the file does not exist.
    """

    if not args.file:
        if sys.stdin.isatty():
            raise ValueError("[id2aant.py] Error: No input provided. Use -h for help.")
        args.file = sys.stdin.read().splitlines()
    else:
        if not os.path.isfile(args.file):
            raise FileNotFoundError(f"[id2aant.py] Error: '{args.file}' does not exist.")

        try:
            args.file = read_ids(args.file)
        except ValueError as e:
            raise e


def make_request_with_retries(url: str, retries=3, delay=1) -> requests.Response:
        """
        Make a request to a given URL with retries in case of failure.

        Parameters:
            url (str): The URL to send the request to.
            retries (int): The maximum number of retries in case of failure (default: 3).
            delay (float): The delay in seconds between retries. Increses exponentially.(default: 1).

        Returns:
            requests.Response: The response object.
        """

        retry_strategy = Retry(
                            total=retries,
                            status_forcelist=[429,500,502,503,504],
                            allowed_methods=["GET", "POST"],
                            backoff_factor=delay,
                            )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        http = requests.Session()
        http.mount("https://", adapter)

        response = http.get(url, timeout=10)

        return response


def handle_uniprot_json_entry(json_entry: str) -> dict:
    """
    Parse UniProt JSON entry and try to acces the important data.

    Args:
        json_entry: UniProt JSON entry.

    Returns:
        entry_data: Dictionary containing the following keys:
            uniprot_accession: UniProt accession.
            uniprot_aa_sequence: UniProt amino acid sequence.
            ena_accession: ENA accession.

    Raises:
        KeyError: If the JSON entry does not contain the expected data.
        ValueError: If the JSON entry does not contain the expected data.
    """

    entry_data = json.loads(json_entry)

    if entry_data["entryType"] == "Inactive":
        raise ValueError(f"[id2aant.py] Warning: {entry_data['primaryAccession']} Uniprot entry is inactive."
                        + " This protein will be skipped.")


    uniprot_accession = entry_data["primaryAccession"]


    for dbcrossref in entry_data["uniProtKBCrossReferences"]:
        if dbcrossref["database"] == "EMBL":
            for prop in dbcrossref["properties"]:
                if prop["key"] == "ProteinId":
                    ena_accession = prop["value"].split(".")[0]
                    break


    uniprot_aa_sequence = entry_data["sequence"]["value"]


    # FIXME: This is a temporary fix. I need to find a better way to handle this.
    if len(uniprot_aa_sequence) and len(ena_accession) and len(uniprot_accession):
        return {"uniprot_accession": uniprot_accession,
                "uniprot_aa_sequence": uniprot_aa_sequence,
                "ena_accession": ena_accession}

    raise ValueError(f"[id2aant.py] Warning: {uniprot_accession} Uniprot entry is missing data."
                    + " This protein will be skipped.")


def get_uniprot_entry_data(uniprot_accession: str) -> dict:
    """
    Get UniProt entry data. Request the JSON entry and parse it.
    If the request fails, entry is inactive or there is a problem parsing the
    JSON, return None.

    Args:
        uniprot_accession: UniProt accession.

    Returns:
        entry_data: Dictionary containing the following keys:
            uniprot_accession: UniProt accession.
            uniprot_aa_sequence: UniProt amino acid sequence.
            ena_accession: ENA accession.

    Raises:
        requests.exceptions.RequestException: If the request fails.
        KeyError: If the JSON entry does not contain the expected data.
        ValueError: If the JSON entry does not contain the expected data.
    """

    url = UNIPROT_JSON_ENTRY.format(code=uniprot_accession)

    response = make_request_with_retries(url)

    if response.status_code != 200:
        raise requests.exceptions.RequestException(f"Warning: request to {url} "
                         + "failed with status code {response.status_code}."
                         + f"Reason: {response.reason}. {uniprot_accession} "
                         + "will be skipped.\n")

    try:
        entry_data = handle_uniprot_json_entry(response.text)
    except KeyError as e:
        raise KeyError(str(e) + f" Url used: {url}\n")

    if not entry_data:
        raise ValueError(f"[id2aant.py] Warning: {uniprot_accession} Uniprot entry is missing data."
                        + f" Url used: {url}\n")

    return entry_data


def get_ena_nucleotide_sequence(ena_accession: str) -> str:
    """
    Get ENA nucleotide sequence.

    Args:
        ena_accession: ENA accession.

    Returns:
        ena_nucleotide_sequence: ENA nucleotide sequence.

    Raises:
        Exception: If the request fails.
    """

    url = ENA_NUCLEOTIDE_FASTA_SEQUENCE.format(code=ena_accession)

    response = make_request_with_retries(url)

    if response.status_code != 200:
        raise requests.exceptions.RequestException(f"[id2aant.py] Warning: request to {url} "
                         + f"failed with status code {response.status_code}."
                         + f"Reason: {response.reason}. {ena_accession} "
                         + "will be skipped.\n")

    if response.text.startswith(">"):
        fasta_sequence = response.text.splitlines()
        return {"header": fasta_sequence[0].strip(),
                "sequence": "".join(fasta_sequence[1:]).replace("\n", "")}

    raise ValueError(f"[id2aant.py] Warning: {url} for {ena_accession} does not return a FASTA sequence.")


def main():

    args = setup_argparse().parse_args()
    try:
        validate_args(args)
    except (ValueError, FileNotFoundError) as e:
        eprint(f"{e}") if not args.quiet else None
        sys.exit(1)


    ids = args.file
    proteins = [Protein(input_id=i) for i in ids]


    for prot in proteins:

        try:
            uniprot_entry_data = get_uniprot_entry_data(prot.input_id)
        except (requests.exceptions.RequestException, KeyError, ValueError) as e:
            eprint(f"{e}") if not args.quiet else None
            continue

        prot.uniprot_accession = uniprot_entry_data["uniprot_accession"]
        prot.uniprot_aa_sequence = uniprot_entry_data["uniprot_aa_sequence"]
        prot.ena_accession = uniprot_entry_data["ena_accession"]


        try:
            ena_nucleotide_sequence = get_ena_nucleotide_sequence(prot.ena_accession)
        except (requests.exceptions.RequestException, ValueError) as e:
            eprint(f"{e}") if not args.quiet else None
            continue

        prot.ena_nucleotide_sequence = ena_nucleotide_sequence["sequence"]


    proteins = [p for p in proteins if p.is_valid()]
    if len(proteins) == 0:
        eprint("[id2aant.py] Error: No remaining proteins to process.")
        sys.exit(1)


    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    with open(f"{args.output}/{args.prefix}aa.fasta", "w") as f:
        f.write("\n".join([f"{p.generate_fasta_header()}\n{p.uniprot_aa_sequence}" for p in proteins]))
    with open(f"{args.output}/{args.prefix}nt.fasta", "w") as f:
        f.write("\n".join([f"{p.generate_fasta_header()}\n{p.ena_nucleotide_sequence}" for p in proteins]))

    if not args.quiet:
        eprint(f"[id2aant.py] Info: Done. {len(proteins)} proteins processed.")
    sys.exit(0)


if __name__ == "__main__":
    main()

