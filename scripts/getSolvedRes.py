#!/usr/bin/env python3

"""
Given a PDB code with chain, and a pdb file, this script will return a `.json`
file containing the solved residues of the protein.

Usage: getSolvedRes.py <pdb_file> <pdb_code>

Author: Joan Lluis Pons Ramon
Email:  <joanlluispons@gmail.com>
"""

import argparse
import json
import os
import sys

from Bio.PDB import PDBParser


__version__ = "0.1.0"


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


def eprint(*args, **kwargs):
    """Print to stderr"""
    print(*args, file=sys.stderr, **kwargs)


def setup_argparse() -> argparse.ArgumentParser:
    """
    Customize and setup argparse.

    Args:
        None

    Returns:
        parser: An argparse.ArgumentParser object.
    """

    fmt = lambda prog: CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(
            prog="getSolvedRes",
            formatter_class=fmt,
            usage="getSolvedRes.py <pdb_file> <pdb_code> [options]",
            description="""
Given a PDB code with chain, and a pdb file, this script will return a `.json`
file containing the solved residues of the protein.

Structure of the output file:
{
    "<residue_number>": {
        "resname": "A",
        "insertion_code": " "
    },
    ...
}
            """,
            epilog="""
Example:
    getSolvedRes.py 1a0m.pdb 1a0m

For more information, visit: <https://github.com/jllpons/mabp>
            """,
        )

    parser.add_argument(
        "pdb_file",
        metavar="FILE",
        type=str,
        help="Path to the PDB file.",
    )
    parser.add_argument(
        "pdb_code",
        metavar="PDB_CODE",
        type=str,
        help="PDB code with chain identifier.",
    )
    parser.add_argument(
        "-o", "--output",
        metavar="PATH",
        type=str,
        default=f"{os.getcwd()}",
        help="Specify the output directory. Will be created if it doesn't exist. Default: current working directory.",
    )
    parser.add_argument(
        "-p", "--prefix",
        metavar="STR",
        type=str,
        default="",
        help="Specify a prefix for the output file."
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """
    Validate the arguments.

    Args:
        args (argparse.Namespace): The parsed arguments.

    Returns:
        None
    """

    if not os.path.isfile(args.pdb_file):
        raise FileNotFoundError(f"[getSolvedRes.py] Error: The file '{args.pdb_file}' does not exist.")

    if len(args.pdb_code) < 5:
        raise ValueError(f"[getSolvedRes.py] Error: The PDB code '{args.pdb_code}' is not valid.")

    if not os.path.isdir(args.output):
        os.makedirs(args.output)


def get_solved_residues_from_pdb(pdb_code, pdb_file):
    """
    Retrieves the solved residues from a PDB file.

    Args:
        pdb_code (str): The PDB code of the structure with the chain.
        pdb_file (str): The path to the PDB file.

    Returns:
        dict: A dictionary containing information about the solved residues.
            - keys: The sequence identifier of the residue.
            - values:
                - 'resname': The one-letter code of the residue.
                - 'insertion_code': The insertion code of the residue.
    """


    structure = PDBParser().get_structure(pdb_code, pdb_file)

    # The Structure object follows the SMCRA (Structure/Model/Chain/Residue/Atom) architecture:
    model = structure[0]
    # Chain of the given reference PDB structure.
    chain = model[pdb_code[-1]]
    # Residues (disordered residues not included) are classes that cotain: (id, resname, segid)
    residues = chain.get_residues()

    d3to1 = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
        }

    solved_residues = {}

    for r in residues:
        # If the residue it's actually an aminoacid (could be water, an artefact...):
        if r.resname in d3to1:
            # `id` is a tuple that contains contains:
            #     hetero-flag[0], sequence identifier[1] and insertion code[2].
            # `resname` is a 3 letter code of the name of the residue
            solved_residues[r.id[1]] = {
                "resname": d3to1[r.resname],
                "insertion_code": r.id[2],
            }

    return solved_residues


def main():

        args = setup_argparse().parse_args()

        try:
            validate_args(args)
        except (FileNotFoundError, ValueError) as e:
            eprint(f"{e}")
            sys.exit(1)

        solved_res_dict = get_solved_residues_from_pdb(args.pdb_code, args.pdb_file)

        output_file = f"{args.output}/{args.prefix}{args.pdb_code}_solved_residues.json"

        with open(output_file, "w") as f:
            json.dump(solved_res_dict, f, indent=4)

        eprint(f"[getSolvedRes.py] Info: Process finished. The solved residues have been saved to '{output_file}'.")
        sys.exit(0)


if __name__ == "__main__":
    main()

