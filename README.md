SVTyper
=======
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/hall-lab/svtyper/master/LICENSE)
[![Build Status](https://travis-ci.org/hall-lab/svtyper.svg?branch=master)](https://travis-ci.org/hall-lab/svtyper)

Bayesian genotyper for structural variants

## Overview

SVTyper performs breakpoint genotyping of structural variants (SVs) using whole genome sequencing data. Users must supply a VCF file of sites to genotype (which may be generated by [LUMPY](https://github.com/arq5x/lumpy-sv)) as well as a BAM/CRAM file of Illumina paired-end reads aligned with [BWA-MEM](https://github.com/lh3/bwa). SVTyper assesses discordant and concordant reads from paired-end and split-read alignments to infer genotypes at each site. Algorithm details and benchmarking are described in [Chiang et al., 2015](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3505.html).

![NA12878 heterozygous deletion](etc/het.png?raw=true "NA12878 heterozygous deletion")

## Installation

Requirements:
- Python 2.7.x

### Install via `pip`

    pip install git+https://github.com/hall-lab/svtyper.git

`svtyper` depends on [pysam][0] _(version 0.15.0 or newer)_, [numpy][1], and [scipy][2]; `svtyper-sso` additionally depends on [cytoolz][7]. If the dependencies aren't already available on your system, `pip` will attempt to download and install them.

## `svtyper` vs `svtyper-sso`

`svtyper` is the original implementation of the genotyping algorithm, and works with multiple samples.  `svtyper-sso` is an alternative implementation of `svtyper` that is optimized for genotyping a single sample.  `svtyper-sso` is a parallelized implementation of `svtyper` that takes advantage of multiple CPU cores via the [multiprocessing][8] module. `svtyper-sso` can offer a 2x or more speedup (depending on how many CPU cores used) in genotyping a single sample. **_NOTE: svtyper-sso is not yet stable. There are minor logging differences between the two and svtyper-sso may exit with an error prematurely when processing CRAM files._**

## Example Usage

### `svtyper`

#### As a Command Line Python Script

```bash
svtyper \
    -i sv.vcf \
    -B sample.bam \
    -l sample.bam.json \
    --read_names_out \
    > sv.gt.vcf
```

#### As a Python Library

```python
import svtyper.classic as svt

input_vcf = "/path/to/input.vcf"
input_bam = "/path/to/input.bam"
library_info = "/path/to/library_info.json"
output_vcf = "/path/to/output.vcf"

with open(input_vcf, "r") as inf, open(output_vcf, "w") as outf:
    svt.sv_genotype(bam_string=input_bam,
                    vcf_in=inf,
                    vcf_out=outf,
                    min_aligned=20,
                    split_weight=1,
                    disc_weight=1,
                    num_samp=1000000,
                    lib_info_path=library_info,
                    debug=False,
                    alignment_outpath=None,
                    ref_fasta=None,
                    sum_quals=False,
                    max_reads=None,
                    read_names_out=True,
                    both_sides=False)

# Results will be inside the /path/to/output.vcf file
```

### New Features

- **Python version** works with python 3+, originally forked from https://github.com/brentp/svtyper

- **Read counting** Counts are the number of reads supporting an SV, not rounded mapq values. This change was made so that counts could be outputted for single cell data were typically only 1 read may support an SV. Original implementation of SVtyper output rounded mapq values, so if there was 1 read with SV support it would be rounded to 0.

- **Read Names Output**: Using the `--read_names_out` flag will create a companion file (*.readnames) containing the names of reads supporting each structural variant. The output format is:
  ```
  sv_id,sample,split_reads,span_reads,clip_reads
  SV1,SAMPLE1,read1,read2,read3
  ```

- **Duplicates** The `--keep_duplicates` flag will keep duplicates for read counting (default: False).

- **Both Sides Option**: The `--both_sides` flag requires reads to align to both sides of a breakpoint to be counted as support.

### `svtyper-sso`

#### As a Command Line Python Script

```bash
svtyper-sso \
    --core 2 # number of cpu cores to use \
    --batch_size 1000 # number of SVs to process in a single batch (default: 1000) \
    --max_reads 1000 # skip genotyping if SV contains valid reads greater than this threshold (default: 1000) \
    -i sv.vcf \
    -B sample.bam \
    -l sample.bam.json \
    > sv.gt.vcf
```

#### As a Python Library

```python
import svtyper.singlesample as sso

input_vcf = "/path/to/input.vcf"
input_bam = "/path/to/input.bam"
library_info = "/path/to/library_info.json"
output_vcf = "/path/to/output.vcf"

with open(input_vcf, "r") as inf, open(output_vcf, "w") as outf:
    sso.sso_genotype(bam_string=input_bam,
                     vcf_in=inf,
                     vcf_out=outf,
                     min_aligned=20,
                     split_weight=1,
                     disc_weight=1,
                     num_samp=1000000,
                     lib_info_path=library_info,
                     debug=False,
                     alignment_outpath=None,
                     ref_fasta=None,
                     sum_quals=False,
                     max_reads=1000,
                     cores=2,
                     batch_size=1000)

# Results will be inside the /path/to/output.vcf file
```

## Development

Requirements:
- Python 2.7 or newer
- GNU Make
- [virtualenv][3] _(or [conda][4] for [anaconda][5] or [miniconda][6] users)_

### Setting Up a Development Environment

#### Using `virtualenv`

    git clone https://github.com/hall-lab/svtyper.git
    cd svtyper
    virtualenv myvenv
    source myvenv/bin/activate
    pip install -e .
    <add, edit, or delete code>
    make test

    # when you're finished with development
    git push <remote-name> <branch>
    deactivate
    cd .. && rm -rf svtyper

#### Using `conda`

    git clone https://github.com/hall-lab/svtyper.git
    cd svtyper
    conda create --channel bioconda --name mycenv pysam numpy scipy cytoolz # type 'y' when prompted with "proceed ([y]/n)?"
    source activate mycenv
    pip install -e .
    <add, edit, or delete code>
    make test


    # when you're finished with development
    git push <remote-name> <branch>
    source deactivate
    cd .. && rm -rf svtyper
    conda remove --name mycenv --all

## Troubleshooting

Many common issues are related to abnormal insert size distributions in the BAM file. SVTyper provides methods to assess and visualize the characteristics of sequencing libraries.

Running SVTyper with the `-l` flag creates a JSON file with essential metrics on a BAM file. SVTyper will sample the first N reads for the file (1 million by default) to parse the libraries, read groups, and insert size histograms. This can be done in the absence of a VCF file.
```
svtyper \
    -B my.bam \
    -l my.bam.json
```

The [lib_stats.R](scripts/lib_stats.R) script produces insert size histograms from the JSON file
```
scripts/lib_stats.R my.bam.json my.bam.json.pdf
```
![Insert size histogram](etc/my.bam.json.png?raw=true "Insert size histogram")


## Citation

C Chiang, R M Layer, G G Faust, M R Lindberg, D B Rose, E P Garrison, G T Marth, A R Quinlan, and I M Hall. SpeedSeq: ultra-fast personal genome analysis and interpretation. Nat Meth 12, 966–968 (2015). doi:10.1038/nmeth.3505.

http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3505.html

[0]: https://github.com/pysam-developers/pysam
[1]: http://www.numpy.org/
[2]: https://www.scipy.org/
[3]: https://github.com/pypa/virtualenv
[4]: https://conda.io/docs/index.html
[5]: https://docs.continuum.io/anaconda/
[6]: https://conda.io/miniconda.html
[7]: https://github.com/pytoolz/cytoolz
[8]: https://docs.python.org/2/library/multiprocessing.html
