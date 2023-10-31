import os.path

import pytest
import vcf

from permute_vcf.contigs import ContigsTable, get_contig_position


@pytest.mark.parametrize(
    "offset,raw_position,chrom,chrom_position",
    [
        (0, 50, "chr1", 50),
        (0, 100, "chr1", 100),
        (0, 101, "chr2", 1),
        (0, 150, "chr2", 50),
        (0, 500, "chr3", 200),
        (40, 50, "chr1", 50),
        (150, 50, "chr2", 50),
        (150, 51, "chr3", 1),
    ],
)
def test_get_contig_position(offset, raw_position, chrom, chrom_position):
    contigs_table = make_contigs_table()
    positions, offsets = contigs_table._make_offsets(offset)

    assert get_contig_position(offsets, raw_position) == (chrom, chrom_position)


def test_zero_offset():
    offsets = make_contigs_table().zero_offsets

    assert offsets[0] == ("chr1", 0)
    assert offsets[5] == ("chr6", 1500)
    assert len(offsets) == 10


def test_small_offset():
    offset = 40
    contigs_table = make_contigs_table()
    positions, offsets = contigs_table._make_offsets(offset)

    assert positions == contigs_table.total_length - offset * contigs_table.contig_count
    assert len(offsets) == 10
    assert offsets[5] == ("chr6", sum((i + 1) * 100 - offset for i in range(5)))


def test_large_offset():
    offset = 150
    contigs_table = make_contigs_table()
    positions, offsets = contigs_table._make_offsets(offset)

    assert positions == sum(
        (i + 1) * 100 - offset for i in range(10) if (i + 1) * 100 > offset
    )
    assert len(offsets) == 9
    assert offsets[5] == ("chr7", sum(i * 100 - offset for i in range(2, 7)))


def test_sample(offset):
    contigs_table = make_contigs_table()
    contigs_table.total_length * 100

    contigs_table.sample(
        min_3p_dist=offset, permutations=contigs_table.total_length * 100
    )


def test_contigs_table():
    contigs_table = make_contigs_table()

    assert contigs_table.total_length == sum((i + 1) * 100 for i in range(10))
    assert contigs_table.contig_lengths["chr6"] == 600
    assert contigs_table.contig_count == 10


def make_contigs_table():
    vcf_path = os.path.join(os.path.dirname(__file__), "data", "test.vcf")
    reader = vcf.Reader(filename=vcf_path)
    return ContigsTable(reader)
