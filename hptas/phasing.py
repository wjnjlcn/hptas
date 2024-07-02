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


import os
import itertools

from typing import List, Set, Dict, Tuple, Optional
from datetime import datetime

import numpy as np
import pandas
import pysam

import HTSeq
from HTSeq import GenomicInterval, GenomicArrayOfSets

from .annotation import Segment, SNP

import stan
import arviz as az


KMER_LEN = 32


def log_msg(msg: str):
    now = datetime.now()
    time = '{}-{:02}-{:02} {:02}:{:02}:{:02}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)
    print('[{}] {}'.format(time, msg))


def encode_kmer(kmer: str) -> Optional[int]:
    kmer_char_list = list(kmer)
    bin_char_list = ['0b']
    for char in kmer_char_list:
        bin_char = None
        if char == 'A' or char == 'a':
            bin_char = '00'
        elif char == 'C' or char == 'c':
            bin_char = '01'
        elif char == 'G' or char == 'g':
            bin_char = '10'
        elif char == 'T' or char == 't':
            bin_char = '11'
        if bin_char is None:
            return None
        bin_char_list.append(bin_char)
    encoded_kmer = int(''.join(bin_char_list), 2)
    # print(encoded_kmer, '{0:b}'.format(encoded_kmer))
    return encoded_kmer


class SNPAllele:
    def __init__(self, snp: SNP, number):
        self.snp = snp
        assert number == 0 or number == 1
        self.number = number

    def nt(self) -> str:
        return self.snp.alleles[self.number]

    def generate_code(self) -> int:
        code = bin(self.snp.idx)
        code += str(self.number)
        return int(code, 2)

    def __repr__(self):
        return '{}:{}'.format(self.snp.id, self.nt())


class SNPAlleleInfo:
    def __init__(self, snp_idx: int, snp_id: str, nt_num: int, nt: str):
        self.gene_id = None  # type: Optional[str]
        self.gene_name = None  # type: Optional[str]
        self.snp_idx = snp_idx
        self.snp_id = snp_id
        self.nt_num = nt_num
        self.nt = nt


class PhasingInfo:
    def __init__(self, gene_id: str, gene_name: str,
                 snp_1_idx: int, snp_1_id: str,
                 snp_2_idx: int, snp_2_id: str,
                 type_1_count: int=0, type_2_count: int=0):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.snp_1_idx = snp_1_idx
        self.snp_1_id = snp_1_id
        self.snp_2_idx = snp_2_idx
        self.snp_2_id = snp_2_id
        self.type_1_count = type_1_count
        self.type_2_count = type_2_count


def get_idx_from_snp_allele_code(snp_allele_code: int) -> int:
    if snp_allele_code == 0:
        code = '0b00'
    elif snp_allele_code == 1:
        code = '0b01'

    else:
        code = bin(snp_allele_code)
    code = list(code)
    idx = int(''.join(code[:-1]), 2)
    return idx


class SNPKmerHTItem:
    def __init__(self, snp_allele_code: int, has_other_items: bool):
        self.snp_allele_code = snp_allele_code
        self.has_other_items = has_other_items

    def generate_code(self) -> int:
        code = bin(self.snp_allele_code)
        code += str(1 if self.has_other_items else 0)
        return int(code, 2)

def create_SNPKmerHTItem_from_code(item_code: int) -> SNPKmerHTItem:
    if item_code == 0:
        code = '0b00'
    elif item_code == 1:
        code = '0b01'

    else:
        code = bin(item_code)
    code = list(code)
    snp_allele_code, has_other_items = int(''.join(code[:-1]), 2), False if int(code[-1]) == 0 else True
    return SNPKmerHTItem(snp_allele_code, has_other_items)


class PhasingHT:
    def __init__(self):
        self.snp_kmer_ht = {}  # type: Dict[int, int]
        self.alt_snp_kmer_ht = {}  # type: Dict[int, List[int]]
        self.result_ht = {}  # type: Dict[Tuple[int, int], int]

    def clear_without_result_ht(self):
        self.snp_kmer_ht = {}
        self.alt_snp_kmer_ht = {}

    def add_gene_snps_for_analysis(self, snps: List[SNP]):
        snps = sorted(snps, key=lambda x: x.pos)
        snp_allele_codes_list = [(SNPAllele(snp, 0).generate_code(),
                                  SNPAllele(snp, 1).generate_code()) for snp in snps]
        for i in range(len(snps)-1):
            c_l = snp_allele_codes_list[i]
            c_r = snp_allele_codes_list[i+1]
            for (l, r) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                self.result_ht[(c_l[l], c_r[r])] = 0

    def get_snp_allele_codes_of_kmer_code(self, encoded_kmer: int) -> List[int]:
        if encoded_kmer is None:
            return []

        try:
            snp_kmer_ht_item_code = self.snp_kmer_ht[encoded_kmer]
            snp_kmer_ht_item = create_SNPKmerHTItem_from_code(snp_kmer_ht_item_code)
            snp_allele_codes = [snp_kmer_ht_item.snp_allele_code]
            if snp_kmer_ht_item.has_other_items:
                snp_allele_codes += self.alt_snp_kmer_ht[encoded_kmer]
            return snp_allele_codes

        except KeyError:
            return []

    def update_result(self, snp_allele_code_pair: Tuple[int, int]):
        try:
            self.result_ht[snp_allele_code_pair] += 1
            # print(snp_allele_code_pair, self.result_ht[snp_allele_code_pair])
        except KeyError:
            pass

    def insert_snp_kmer_with_snp_alleles(self, snp_kmer: str, snp_alleles: List[SNPAllele]):
        # print(snp_kmer, len(snp_kmer), snp_alleles)
        encoded_snp_kmer = encode_kmer(snp_kmer)
        if encoded_snp_kmer is None:
            return

        snp_allele_codes = [snp_allele.generate_code() for snp_allele in snp_alleles]
        snp_allele_codes = list(set(snp_allele_codes))  # Remove duplicates.
        # print(snp_allele_codes)
        try:
            snp_kmer_ht_item_code = self.snp_kmer_ht[encoded_snp_kmer]
            snp_kmer_ht_item = create_SNPKmerHTItem_from_code(snp_kmer_ht_item_code)
            other_snp_allele_codes = [c for c in snp_allele_codes if c != snp_kmer_ht_item.snp_allele_code]
        except KeyError:
            snp_kmer_ht_item = SNPKmerHTItem(snp_allele_codes[0], False)
            self.snp_kmer_ht[encoded_snp_kmer] = snp_kmer_ht_item.generate_code()
            other_snp_allele_codes = snp_allele_codes[1:]
        # print(snp_allele_codes, other_snp_allele_codes)
        if len(other_snp_allele_codes) > 0:
            snp_kmer_ht_item.has_other_items = True
            self.snp_kmer_ht[encoded_snp_kmer] = snp_kmer_ht_item.generate_code()
            self._insert_encoded_snp_kmer_with_other_snp_allele_codes(encoded_snp_kmer, other_snp_allele_codes)

    def _insert_encoded_snp_kmer_with_other_snp_allele_codes(self, encoded_snp_kmer: int,
                                                             other_snp_allele_codes: List[int]):
        try:
            alt_snp_allele_codes = self.alt_snp_kmer_ht[encoded_snp_kmer]
        except KeyError:
            alt_snp_allele_codes = []
            self.alt_snp_kmer_ht[encoded_snp_kmer] = alt_snp_allele_codes
        for other_snp_allele_code in other_snp_allele_codes:
            if other_snp_allele_code not in alt_snp_allele_codes:
                alt_snp_allele_codes.append(other_snp_allele_code)


SNP_FILENAME = 'snps.txt'
PHASING_COUNTS_FILENAME = 'counts.txt'
PHASING_RESULT_FILENAME = 'result.csv'

PHASING_CODE = """
    data {
      int<lower=0> N;
      int<lower=0> H;
    }
    parameters {
      real<lower=0, upper=1> P;
    }
    transformed parameters {
    }
    model {
      P ~ beta(0.5, 0.5);
      H ~ binomial(N, P);
    }
"""

class PhasingError(Exception):
    def __init__(self, msg):
        super().__init__(msg)


class PhasingPathError(PhasingError):
    """Path error."""
    def __init__(self, path, msg):
        super().__init__('[PATH: {}] {}'.format(path, msg))


class PhasingCount:
    def __init__(self, ref_fa_path: str, result_dir_path: str, error_kmer_path: Optional[str]=None):
        self.pht = PhasingHT()

        self.ref_fa_path = ref_fa_path
        self.result_dir_path = result_dir_path
        self.error_kmer_codes = set()

        self.snps_count = 0

        if os.path.exists(self.snp_file_path):
            raise PhasingPathError(self.snp_file_path, 'SNP file already exists in the output directory.')

        self._init_error_kmers(error_kmer_path)

    @property
    def snp_file_path(self) -> str:
        return os.path.join(self.result_dir_path, SNP_FILENAME)

    @property
    def counts_file_path(self) -> str:
        return os.path.join(self.result_dir_path, PHASING_COUNTS_FILENAME)

    def _init_error_kmers(self, error_kmer_path: Optional[str]):
        if error_kmer_path is not None:
            if not os.path.exists(error_kmer_path):
                raise PhasingPathError(error_kmer_path, 'The error k-mer file does not exist.')

            with open(error_kmer_path) as f:
                for line in f:
                    kmer = line.strip('\n')
                    if len(kmer) != KMER_LEN:
                        continue
                    kmer_code = encode_kmer(kmer)
                    self.error_kmer_codes.add(kmer_code)

    def add_gene_with_snps_for_analysis(self, gene: Segment, snps: List[SNP]):
        # Sort SNPs by position.
        snps = sorted(snps, key=lambda x: x.pos)
        # Create SNP idx for each SNP.
        for i, snp in enumerate(snps):
            snp.idx = i + self.snps_count
            snp.gene_id = gene.id
            snp.gene_name = gene.gene_name
        self.snps_count += len(snps)

        with open(self.snp_file_path, 'a') as snp_file:
            for snp in snps:
                snp_file.write('\t'.join([snp.gene_id, snp.gene_name,
                                          str(snp.idx), snp.id, snp.chrom, str(snp.pos),
                                          snp.alleles[0], snp.alleles[1]]) + '\n')

        self.pht.add_gene_snps_for_analysis(snps)

        with pysam.FastaFile(self.ref_fa_path) as ref_fa:
            gene_seq = ref_fa.fetch(reference=gene.iv.chrom, start=gene.iv.start, end=gene.iv.end)

        for iso_num in range(gene.isoforms_count):
            self._extract_snp_kmers_for_isoform_num(iso_num, gene, gene_seq, snps)

    def _finish_and_write_results(self):
        self.pht.clear_without_result_ht()

        snp_idx_map = {}
        with open(self.snp_file_path) as snp_file:
            for line in snp_file:
                cols = line.strip('\n').split('\t')
                snp_gene_id, snp_gene_name, snp_idx, snp_id, snp_chrom, pos, a1, a2 = cols
                snp = SNP(snp_chrom, int(pos), [a1, a2], False)
                snp.gene_id = snp_gene_id
                snp.gene_name = snp_gene_name
                snp.idx = int(snp_idx)
                snp_idx_map[snp.idx] = snp
        with open(self.counts_file_path, 'w') as out:
            out.write('#Gene-ID\tGene-NAME\t')
            for i in range(2):
                out.write('SNP#{}-IDX\tSNP#{}-ID\tSNP#{}-NTNUM\tSNP#{}-NT\t'.format(i+1, i+1, i+1, i+1))
            out.write('Count\n')

            for snp_allele_code_pair, count in self.pht.result_ht.items():
                snp_allele_pair = [self._create_SNPAllele_from_code(c, snp_idx_map) for c in snp_allele_code_pair]
                out.write('{}\t{}\t'.format(snp_allele_pair[0].snp.gene_id, snp_allele_pair[0].snp.gene_name))
                for s in snp_allele_pair:
                    out.write('\t'.join([str(s.snp.idx), s.snp.id, str(s.number), s.nt()]) + '\t')
                out.write('{}\n'.format(count))


    def _create_SNPAllele_from_code(self, snp_allele_code: int, snp_idx_map: Dict[int, SNP]) -> SNPAllele:
        idx = get_idx_from_snp_allele_code(snp_allele_code)
        snp = snp_idx_map[idx]
        allele_num = int(list(bin(snp_allele_code))[-1])
        return SNPAllele(snp, allele_num)

    def _extract_snp_kmers_for_isoform_num(self, n, gene: Segment, gene_seq: str, snps: List[SNP]):
        isoform = gene.isoforms[n]
        ga = GenomicArrayOfSets('auto', False)
        for exon in isoform.exons:
            ga[exon] += 'e'
        for snp in snps:
            ga[snp.iv] += snp

        kmer_window_count = isoform.length - KMER_LEN + 1
        if kmer_window_count <= 0:
            return

        n = 0
        for iv, v in ga.steps():
            if 'e' in v:
                for start in range(iv.start, iv.end):
                    self._extract_snps_kmer_start_from(start, ga, isoform, gene_seq, gene.iv.start)
                    n += 1
                    if n >= kmer_window_count:
                        return

    def _extract_snps_kmer_start_from(self, start: int, ga: GenomicArrayOfSets, isoform, seq: str, seq_start: int):
        iv_range = GenomicInterval(isoform.iv.chrom, start, isoform.iv.end, '.')
        ivs = []
        snps_with_idx = []

        length = 0
        for iv, v in ga[iv_range].steps():
            # print(iv, v)
            if 'e' in v:
                end_flag = False
                e = iv.end
                if length + iv.length >= KMER_LEN:
                    end_flag = True
                    e = iv.start + (KMER_LEN - length)
                selected_iv = GenomicInterval(iv.chrom, iv.start, e, '.')
                length += selected_iv.length
                ivs.append(selected_iv)
                if len(v) > 1:
                    snp = None
                    for elem in v:
                        if isinstance(elem, SNP):
                            snp = elem
                    if snp is not None:
                        snps_with_idx.append((snp, length - 1))
                if end_flag:
                    break
        # print(ivs, snps_with_idx)
        if len(snps_with_idx) > 0:
            selected_seq = ''
            for iv in ivs:
                selected_seq += seq[iv.start - seq_start: iv.end - seq_start]
            mod_array = list(itertools.product([0, 1], repeat=len(snps_with_idx)))
            for mod in mod_array:
                snp_kmer = list(selected_seq)
                snp_alleles = []
                for k, (snp, idx) in enumerate(snps_with_idx):
                    snp_allele = SNPAllele(snp, mod[k])
                    snp_alleles.append(snp_allele)
                    snp_kmer[idx] = snp_allele.nt()

                snp_kmer = ''.join(snp_kmer)
                self.pht.insert_snp_kmer_with_snp_alleles(snp_kmer, snp_alleles)

    def analyze_seqs(self, fq_paths: List[str], paired_end: bool):
        if paired_end:
            assert len(fq_paths) == 2
        # else:
        #     assert len(fq_paths) == 1

        analyzed_frags_report_interval = 100000

        if paired_end:
            log_msg('Analyzing PE Fastq files: {}'.format(fq_paths))
            n = 0
            with HTSeq.FastqReader(fq_paths[0]) as fq_1, HTSeq.FastqReader(fq_paths[1]) as fq_2:
                for seq1, seq2 in zip(fq_1, fq_2):
                    # print(seq1, seq2)
                    self._analyze_paired_end_fq_seqs(seq1, seq2)
                    n += 1
                    if n % analyzed_frags_report_interval == 0:
                        log_msg('{} read pairs analyzed.'.format(n))
        else:
            for fq_path in fq_paths:
                log_msg('Analyzing Fastq file: {}'.format(fq_path))
                n = 0
                with HTSeq.FastqReader(fq_path) as fq:
                    for seq in fq:
                        self._analyze_single_fq_seq(seq)
                        n += 1

                        # For testing.
                        # if n >= 10000:
                        #     break

                        if n % analyzed_frags_report_interval == 0:
                            log_msg('{} reads analyzed.'.format(n))

        self._finish_and_write_results()

    def _analyze_paired_end_fq_seqs(self, fq_seq_1: HTSeq.SequenceWithQualities, fq_seq_2: HTSeq.SequenceWithQualities):
        snp_allele_codes = self._get_snp_allele_codes_of_fq_seq(fq_seq_1)
        for c in self._get_snp_allele_codes_of_fq_seq(fq_seq_2):
            snp_allele_codes.add(c)
        self._update_fragment_snp_allele_codes(snp_allele_codes)

    def _analyze_single_fq_seq(self, fq_seq: HTSeq.SequenceWithQualities):
        snp_allele_codes = self._get_snp_allele_codes_of_fq_seq(fq_seq)
        self._update_fragment_snp_allele_codes(snp_allele_codes)

    def _update_fragment_snp_allele_codes(self, snp_allele_codes: Set[int]):
        if len(snp_allele_codes) < 2:
            return
        snp_allele_codes = sorted(list(snp_allele_codes), key=get_idx_from_snp_allele_code)
        # print(snp_allele_codes)
        for i in range(len(snp_allele_codes)-1):
            snp_allele_code_pair = (snp_allele_codes[i], snp_allele_codes[i+1])
            self.pht.update_result(snp_allele_code_pair)

    def _get_snp_allele_codes_of_fq_seq(self, fq_seq: HTSeq.SequenceWithQualities) -> Set[int]:
        # Fastq read or its Reverse Complement.
        snp_allele_codes = self._get_snp_allele_codes_of_seq(str(fq_seq.seq, 'ascii'))
        # print('+', snp_allele_codes)
        if len(snp_allele_codes) == 0:
            snp_allele_codes = self._get_snp_allele_codes_of_seq(str(fq_seq.get_reverse_complement().seq, 'ascii'))
            # print('-', snp_allele_codes)
        return snp_allele_codes

    def _get_snp_allele_codes_of_seq(self, seq: str) -> Set[int]:
        if len(seq) < KMER_LEN:
            return set()

        n = int(len(seq) / KMER_LEN)

        extract_additional_tail_kmer_flag = False
        if len(seq) % KMER_LEN != 0:
            extract_additional_tail_kmer_flag = True

        kmers = []
        for i in range(n):
            kmer = seq[i*KMER_LEN: (i+1)*KMER_LEN]
            kmers.append(kmer)
        if extract_additional_tail_kmer_flag:
            additional_tail_kmer = seq[len(seq)-KMER_LEN:]
            kmers.append(additional_tail_kmer)

        snp_allele_codes = set()
        for kmer in kmers:
            encoded_kmer = encode_kmer(kmer)
            if encoded_kmer in self.error_kmer_codes:
                continue

            for c in self.pht.get_snp_allele_codes_of_kmer_code(encoded_kmer):
                snp_allele_codes.add(c)

        return snp_allele_codes


class PhasingInference:
    def __init__(self, result_dir_path: str, counts_file_path: Optional[str]=None, result_file_path: Optional[str]=None):
        self.result_dir_path = result_dir_path
        self._counts_file_path = counts_file_path
        self._result_file_path = result_file_path

    @property
    def counts_file_path(self) -> str:
        if self._counts_file_path is not None:
            return self._counts_file_path
        return os.path.join(self.result_dir_path, PHASING_COUNTS_FILENAME)

    @property
    def result_file_path(self) -> str:
        if self._result_file_path is not None:
            return self._result_file_path
        return os.path.join(self.result_dir_path, PHASING_RESULT_FILENAME)

    def bayesian_inference(self):
        counts_d = self.load_counts_data()

        cols = ['Gene-ID', 'Gene-NAME']
        for i in range(2):
            cols += ['SNP#%d-IDX' % (i+1), 'SNP#%d-ID' % (i+1)]
        cols += ['Type#1(00/11)-COUNT', 'Type#2(01/10)-COUNT']
        cols += ['Mean', 'STD', '95%HDI-Width', 'HDI-L', 'HDI-U']

        d = [[] for _ in range(len(cols))]

        for _, phasing_info in counts_d.items():
            mean, std, hdi_l, hdi_u = self._bayesian_inference_for_one(phasing_info.type_1_count, phasing_info.type_2_count)
            for i, v in enumerate([phasing_info.gene_id, phasing_info.gene_name,
                                   phasing_info.snp_1_idx, phasing_info.snp_1_id,
                                   phasing_info.snp_2_idx, phasing_info.snp_2_id,
                                   phasing_info.type_1_count, phasing_info.type_2_count,
                                   mean, std, hdi_u-hdi_l, hdi_l, hdi_u]):
                d[i].append(v)
        df = pandas.DataFrame({col: d[i] for i, col in enumerate(cols)})[cols]
        df.to_csv(self.result_file_path, index=False)

    def load_counts_data(self) -> Dict[Tuple[int, int], PhasingInfo]:
        d = {}  # type: Dict[Tuple[int, int], PhasingInfo]

        with open(self.counts_file_path) as counts_file:
            for line in counts_file:
                if line.startswith('#'):
                    continue
                cols = line.strip('\n').split('\t')
                gene_id = cols[0]
                gene_name = cols[1]
                snp_allele_infos = []
                for i in range(2):
                    snp_idx, snp_id, nt_num, nt = cols[i*4+2:i*4+2+4]
                    snp_allele_info = SNPAlleleInfo(int(snp_idx), snp_id, int(nt_num), nt)
                    snp_allele_info.gene_id = gene_id
                    snp_allele_info.gene_name = gene_name
                    snp_allele_infos.append(snp_allele_info)
                count = int(cols[-1])

                type_1_count = 0
                type_2_count = 0
                nt_num_pair = (snp_allele_infos[0].nt_num, snp_allele_infos[1].nt_num)
                if nt_num_pair in [(0, 0), (1, 1)]:
                    type_1_count = count
                else:
                    type_2_count = count

                key = (snp_allele_infos[0].snp_idx, snp_allele_infos[1].snp_idx)
                try:
                    phasing_info = d[key]
                except KeyError:
                    phasing_info = PhasingInfo(snp_allele_infos[0].gene_id,
                                               snp_allele_infos[0].gene_name,
                                               snp_allele_infos[0].snp_idx,
                                               snp_allele_infos[0].snp_id,
                                               snp_allele_infos[1].snp_idx,
                                               snp_allele_infos[1].snp_id)
                    d[key] = phasing_info

                phasing_info.type_1_count += type_1_count
                phasing_info.type_2_count += type_2_count
        return d


    def _bayesian_inference_for_one(self, type_1_count, type_2_count):
        if type_1_count + type_2_count < 2:
            return np.nan, np.nan, np.nan, np.nan

        phasing_data = {'N': type_1_count + type_2_count, 'H': type_1_count}

        posterior = stan.build(PHASING_CODE, data=phasing_data)
        fit = posterior.sample(num_chains=4, num_samples=1000)
        samples = fit['P'][0]

        mean = np.mean(samples)
        std = np.std(samples)

        hdi_l, hdi_u = az.hdi(samples, hdi_prob=0.95)

        return mean, std, hdi_l, hdi_u
