#!/usr/bin/env python3

"""
From a list of UniProt IDs, this script fetches the amino acid sequences from
UniProt and the nucleotide sequences from ENA. Then, it writes them in two
fasta files. FASTA headers will appear as: `>UniprotAccession_ENAAccession`.
Accepts IDs list from sdin or from file if -f is used.

Usage:  ids2aant.py [-f] < UniProtIDsList.txt > output.fasta

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


__version__ = "0.1.0"


def eprint(*args, **kwargs):
    """Print to stderr"""
    print(*args, file=sys.stderr, **kwargs)


# <https://stackoverflow.com/questions/18275023/dont-show-long-options-twice-in-print-help-from-argparse>
#
# `argparse.RawTextHelpFormatter` is also added to the inheritance
# because I want to be able to use "\n" in the description and epilog.
class CustomHelpFormatter(argparse.RawTextHelpFormatter, argparse.HelpFormatter):
    def __init__(self, prog):
        # Initialize with super from RawTextHelpFormatter and also set our custom widths
        super().__init__(prog, max_help_position=40)

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


def setup_argparse() -> argparse.ArgumentParser:
    """
    Setup argparse
    """

    fmt = lambda prog: CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(
            prog="ids2aant",
            formatter_class=fmt,
            usage="ids2aant.py <FILE> [options]",
            description="""
Given a list of UniProt IDs, fetches the amino acid sequences from UniProt and
the nucleotide sequences from ENA. Then, it writes them in two fasta files.
FASTA headers will appear as: `>UniprotAccession_ENAAccession`.
            """,
            epilog="""
Examples:
    ids2aant.py list.txt -o results
    cat list.txt | ids2aant.py

For more information, visit <https://github.com/jllpons/mabp>.
            """,
            )

    parser.add_argument(
            "file",
            metavar="FILE",
            nargs="?",
            type=str,
            help="File containing UniProt IDs separated by newlines. If not provided, it will be read from stdin.",
            )
    parser.add_argument(
            "-o",
            "--output",
            metavar="DIR",
            type=str,
            default=f"{os.getcwd()}",
            help="Output directory. Default: current directory.",
            )
    parser.add_argument(
            "-V",
            "--version",
            action="version",
            version=f"%(prog)s {__version__}",
            )
    parser.add_argument(
            "-q",
            "--quiet",
            action="store_true",
            help="Do not print anything to stderr.",
            )

    return parser

def main():

    args = setup_argparse().parse_args()

if __name__ == "__main__":
    main()

