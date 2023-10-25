import argparse

import vcf


class ContigsTable:
    def __init__(self, vcf_reader: vcf.Reader):
        self.contig_lengths = {c.id: c.length for c in vcf_reader.contigs.values()}
        self.zero_offsets = self._make_offsets()
        self.total_length = sum(self.contig_lengths.values())
        self.contig_count = len(self.contig_lengths)

    def sample(self, offset: int = 0):
        if offset == 0:
            offsets = self.zero_offsets
            possible_positions = self.total_length
        else:
            offsets = self._make_offsets(offset)
            possible_positions = self.total_length - offset * self.contig_count

        # TODO actually do the sampling

    def _make_offsets(self, size=0) -> list[tuple[str, int]]:
        start_position = 1
        offsets = []
        for contig, length in self.contig_lengths.items():
            offsets.append((contig, start_position))
            start_position += length - size

        return offsets


def parse_args():
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
    contigs_table.sample()


if __name__ == "__main__":
    main()
