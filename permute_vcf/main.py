import argparse
import random

import vcf


class ContigNotFoundError(Exception):
    pass


def get_contig_position(
    offsets: list[tuple[str, int]], position: int
) -> tuple[str, int]:
    """Find the contig coordinates of a given raw position

    When randomly choosing a position in the genome, we first count the
    number of possible positions for this variant in the genome and
    then choose a random integer in the range [1, possible_positions].
    This function translates this single integer into an actual genomic
    coordinate based on the list of offsets.

    Args:
        offsets: a list of (contig, offset) pairs ordered by offset in
            ascending order
        position: raw position to convert to a contig coordinate

    Returns: a tuple of the contig name and 1-based position on that
        contig of the raw input position
    """

    # TODO do a binary search instead of a linear search
    # this will speed things up!
    for contig, offset in offsets:
        if offset < position:
            return (contig, offset)

    raise ContigNotFoundError()  # TODO actually make a nice error message


class ContigsTable:
    """table of contig lengths with convenience functions

    This class holds a dictionary of contig lengths, and also allows
    randomly sampling positions from the genome based on this
    information.
    """

    def __init__(self, vcf_reader: vcf.Reader):
        """Create a ContigsTable based on vcf header

        Args:
            vcf_reader: a Reader corresponding to a vcf file whose
                header contains the contig lengths to be stored in
                this object
        """
        self.contig_lengths = {c.id: c.length for c in vcf_reader.contigs.values()}
        self.total_length, self.zero_offsets = self._make_offsets()
        self.contig_count = len(self.contig_lengths)

    def sample(self, min_3p_dist: int = 0) -> tuple[str, int]:
        """Choose a random position in the genome

        This function randomly chooses a valid coordinate position in
        the genome using a uniform distribution.

        Args:
            min_3p_dist: the distance from the end of contigs that the
                chosen position must be. For example, if we are trying
                to randomly choose a start position position for a 2kb
                deletion, it cannot start less than 2kb from the 3' end
                of the chromosome, because a VCF always uses the
                position of the 5' end of a variant as the official
                start position of that variant.

        Returns: a random position in the genome, represented as a
            tuple (contig name, 1-based position on contig)
        """
        if min_3p_dist == 0:
            offsets = self.zero_offsets
            possible_positions = self.total_length
        else:
            possible_positions, offsets = self._make_offsets(min_3p_dist)

        position = random.randrange(1, possible_positions + 1)
        return get_contig_position(offsets, position)

    def _make_offsets(self, min_3p_dist=0) -> tuple[int, list[tuple[str, int]]]:
        """create an offsets table for sampling

        In order to randomly choose a position in the genome, two
        things are necessary:
            1. the number of possibile positions in the genome, so
               that we can randomly sample from the integer range
               [1, possible positions]
            2. a way to convert a number sampled from above range to
               an actual genomic coordinate

        This function provides both.

        Args:
            min_3p_dist: see description of this argument in sample()

        Returns:
            A tuple (possible_positions, offsets_list).
            possible_positions is the total number of possible starting
            positions for this variant in the genome. offsets_list can
            be used by get_contig_position to convert an integer in the
            range [1, possible_positions] into a genomic coordinate.
        """
        start_position = 1
        offsets = []
        for contig, length in self.contig_lengths.items():
            if min_3p_dist <= length:  # can't have an SV bigger than the contig!
                offsets.append((contig, start_position))
                start_position += length - min_3p_dist

        return (start_position - 1, offsets)


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

    contigs_table = ContigsTable(args.input_vcf)
    print(contigs_table.sample())

    # TODO still need to actually loop through the VCF and permute it etc.


if __name__ == "__main__":
    main()
