from os.path import dirname, join
from unittest.mock import patch

import pytest
import vcf

from permute_vcf.main import main


def test_main_help(capsys):
    with patch("sys.argv", ["permute-vcf", "--help"]):
        with pytest.raises(SystemExit):
            main()

    out, _ = capsys.readouterr()
    assert "Given a VCF file" in out


PERMUTATIONS = 10


def test_main_simple(capsys, tmpdir):
    with patch(
        "sys.argv",
        [
            "permute-vcf",
            "-n",
            str(PERMUTATIONS),
            "-o",
            str(tmpdir),
            join(dirname(__file__), "data", "test.vcf"),
        ],
    ):
        main()

    for i in range(PERMUTATIONS):
        reader = vcf.Reader(filename=join(tmpdir, f"{i}.vcf.gz"))
        for record in reader:
            assert record.POS <= reader.contigs[record.CHROM].length
