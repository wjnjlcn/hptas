# This file is part of HPTAS.
#
# HPTAS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HPTAS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HPTAS.  If not, see <https://www.gnu.org/licenses/>.
#
# Author: Jianan Wang
#


import sys
import getopt

from typing import Optional

from .annotation import generate_annotation, AnnotationDB
from .phasing import PhasingCount, PhasingInference


USAGE = '''
Usage:  hptas <command> [options]

Commands:
    anno        build annotation database
    phasing     compute phasing counts and bayesian inference results
'''

ANNO_USAGE = '''
Usage:  hptas anno [options] OUT_DIR
Options:
    -g, --gff FILE          transcriptome annotation GFF3 file
    -v, --vcf FILE          SNP VCF file
    -s, --sample STR        sample name in the VCF file
    -p, --add-chr-prefix    add "chr" prefix to the chrom column of VCF file
'''


PHASING_USAGE = '''
Usage:  hptas phasing [options] READ_FILE...
Options:
    -f, --ref FILE      reference genome fasta file
    -a, --anno-dir DIR  directory of annotation
    -o, --out-dir DIR   output dir
    -p, --pe            input paired-end files
    -e, --error FILE    error k-mer file
'''


def usage(command_idx: Optional[int] = None):
    if command_idx is None:
        print(USAGE)
    else:
        print([ANNO_USAGE, PHASING_USAGE][command_idx])


def main():
    argvs = sys.argv[1:]
    if len(argvs) == 0:
        usage()
        sys.exit(2)

    if argvs[0] == 'anno':
        try:
            opts, args = getopt.gnu_getopt(argvs[1:], 'g:v:s:ph',
                                           ['gff=', 'vcf=', 'sample=', 'add-chr-prefix', 'help'])
        except getopt.GetoptError as error:
            print(error)
            usage(0)
            sys.exit(2)

        gff_path = None
        vcf_path = None
        vcf_sample = None
        add_chrom_prefix = False

        for o, a in opts:
            if o in ('-h', '--help'):
                usage(0)
                sys.exit()
            if o in ('-g', '--gff'):
                gff_path = a
            elif o in ('-v', '--vcf'):
                vcf_path = a
            elif o in ('-s', '--sample'):
                vcf_sample = a
            elif o in ('-p', '--add-chr-prefix'):
                add_chrom_prefix = True

        for v in (gff_path, vcf_path, vcf_sample):
            if v is None:
                usage(0)
                sys.exit(2)

        if len(args) != 1:
            usage(0)
            sys.exit(2)

        out_dir = args[0]

        generate_annotation(gff_path, vcf_path, vcf_sample, add_chrom_prefix, out_dir)

    elif argvs[0] == 'phasing':
        try:
            opts, args = getopt.gnu_getopt(argvs[1:], 'f:a:o:pe:h',
                                           ['ref=', 'anno-dir=', 'out-dir=', 'pe', 'error=', 'help'])
        except getopt.GetoptError as error:
            print(error)
            usage(1)
            sys.exit(2)

        fa_path = None
        anno_dir = None
        out_dir = None
        paired_end = False
        error_kmer_path = None

        for o, a in opts:
            if o in ('-h', '--help'):
                usage(1)
                sys.exit()
            if o in ('-f', '--ref'):
                fa_path = a
            elif o in ('-a', '--anno-dir'):
                anno_dir = a
            elif o in ('-o', '--out-dir'):
                out_dir = a
            elif o in ('-p', '--pe'):
                paired_end = True
            elif o in ('-e', '--error'):
                error_kmer_path = a

        for v in (fa_path, anno_dir, out_dir):
            if v is None:
                usage(1)
                sys.exit(2)

        if len(args) < 1:
            usage(1)
            sys.exit(2)

        read_files = args

        if paired_end and len(args) != 2:
            print('Number of paired-end read files should be 2.')
            usage(1)
            sys.exit(2)

        phasing = PhasingCount(fa_path, out_dir, error_kmer_path=error_kmer_path)

        with AnnotationDB(anno_dir) as anno_db:
            task_ids = anno_db.get_all_task_ids()
            for task_id in task_ids:
                task = anno_db.get_task(task_id)
                # print(task.id, task.segment.gene_name, len(task.snps))
                phasing.add_gene_with_snps_for_analysis(task.segment, task.snps)

        phasing.analyze_seqs(read_files, paired_end)

        phasing = PhasingInference(out_dir)
        phasing.bayesian_inference()

    else:
        usage()
        sys.exit(2)
