# HPTAS

HPTAS is a library for haplotype phasing using transcriptome
information and RNA-seq data.

## Installation

HPTAS is released as a Python package which requires **Python 3.7** 
or higher to be installed on the computer.

To install HPTAS, first use `pip` to install dependent packages:
```shell
pip install --user numpy pandas
pip install --user pysam htseq
pip install --user stan arviz
```

then, use `pip` to install HPTAS:
```shell
pip install --user hptas
```

# Example

To demonstrate the use of HPTAS, an example is provided, it can be downloaded at: [example_data.zip](https://github.com/wjnjlcn/hptas/raw/main/example_data.zip).

Unzip the file and enter the directory `example_data`.

**Build annotation database**
```shell
hptas anno -g test.hg19.gff3 -v test.snp.vcf -s NA12878 -p anno
```

**Haplotype phasing**
```shell
hptas phasing -f chr1.fa -a anno -o results -p -e error.kmers.txt test.1.fq.gz test.2.fq.gz
```