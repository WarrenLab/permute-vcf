"""permute a VCF file.

Given a VCF file, create a series of new VCF files containing the same
variants, but in random positions.
"""

import argparse
import gzip
from os.path import join

import vcf

from permute_vcf.contigs import ContigsTable


def parse_args():
    """Parse command-line arguments"""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "input_vcf",
        type=lambda f: vcf.Reader(filename=f),
        help="input vcf(.gz) file to permute",
    )
    parser.add_argument(
        "-n",
        "--permutations",
        type=int,
        help="number of permutations to run",
        default=1,
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        help="directory where permutations should be output",
        default="permutations/",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # prepare a contigs table for sampling by reading the header
    contigs_table = ContigsTable(args.input_vcf)

    # open a VCF writer for each permuation
    out_vcfs = []
    for permutation in range(args.permutations):
        out_vcfs.append(
            vcf.Writer(
                gzip.open(join(args.output_directory, f"{permutation}.vcf.gz"), "wt"),
                args.input_vcf,  # this copies input header to output files
            )
        )

    # actually loop through the VCF and permute it etc.
    for record in args.input_vcf:
        sv_size = min(len(alt) for alt in record.ALT) - len(record.REF)
        if sv_size < 0:  # deletion
            min_3p_dist = -sv_size
        else:
            min_3p_dist = 0

        random_positions = contigs_table.sample(
            min_3p_dist=min_3p_dist, permutations=args.permutations
        )

        for permutation, random_position in enumerate(random_positions):
            record.CHROM, record.POS = random_position
            out_vcfs[permutation].write_record(record)


if __name__ == "__main__":
    main()
