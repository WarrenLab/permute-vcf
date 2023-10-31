# permute-vcf

randomly permute the locations of variants in a vcf

In order to perform a permutations test on a set of variants, you first need to
permute the locations of the variants a bunch of times. This program does just
that: it takes as input a VCF file containing the actual locations of variants
and outputs a set of randomized versions of the VCF file, each containing the
same variants but in random positions.

## Installation

Due to a [bug][pyvcf-bug] in [pyvcf][pyvcf], this package only supports python
3.9 and 3.10 at the moment. That bug should be fixed in the next release of
pyvcf3, but in the meantime, you'll need to make a python3.9 environment. So to
install:

```bash
conda create -n permute-vcf python=3.9
conda activate permute-vcf
git clone https://github.com/WarrenLab/permute-vcf.git
cd permute-vcf
pip install .
```

[pyvcf-bug]: https://github.com/dridk/PyVCF3/issues/6 [pyvcf]:
https://github.com/dridk/PyVCF3

## Usage

To create 10 permutations of a vcf file `test.vcf.gz` and write them to the
`permutations/` subdirectory, you can run a command like this:

```bash
permute-vcf -n 10 -o permutations/ test.vcf.gz
```

After that, your output directory will have files `0.vcf.gz` through
`9.vcf.gz`, each containing a different permutation of the input data.

The help dialog for the program has more details about parameters etc.
