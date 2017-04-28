# coding=utf-8
# after we get RNA-seq bed file, we want to know the gene expression value.
# 下面是Refseq gene数据库记录的基因的转录谱id的格式：
# -------------------------------------------------------------------------------------------
# field	example	SQL type	info	description
# bin	637	smallint(5) unsigned	range	Indexing field to speed chromosome range queries.
# name	NM_021010	varchar(255)	values	Name of gene (usually transcript_id from GTF)
# chrom	chr8	varchar(255)	values	Reference sequence chromosome or scaffold
# strand	-	char(1)	values	+ or - for strand
# txStart	6912828	int(10) unsigned	range	Transcription start position
# txEnd	6914259	int(10) unsigned	range	Transcription end position
# cdsStart	6912952	int(10) unsigned	range	Coding region start
# cdsEnd	6914219	int(10) unsigned	range	Coding region end
# exonCount	2	int(10) unsigned	range	Number of exons
# exonStarts	6912828,6914047,	longblob	 	Exon start positions
# exonEnds	6913065,6914259,	longblob	 	Exon end positions
# score	0	int(11)	range	score
# name2	DEFA5	varchar(255)	values	Alternate name (e.g. gene_id from GTF)
# cdsStartStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
# cdsEndStat	cmpl	enum('none', 'unk', 'incmpl', 'cmpl')	values	enum('none','unk','incmpl','cmpl')
# exonFrames	1,0,	longblob	 	Exon frame {0,1,2}, or -1 if no frame for exon
# --------------------------------------------------------------------------------------------------------
from bisect import bisect_left, bisect_right
import os
import numpy as np
from time import time


class RefGene(object):
    """
    bin, name, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_count, exon_starts, exon_ends, score, name2,
    cds_start_stat, cds_end_stat, exon_frames
    """
    def __init__(self):
        self.bin, self.name, self.chrom, self.strand = None, '', '', ''
        self.tx_start, self.tx_end, self.cds_start, self.cds_end = None, None, None, None
        self.exon_count, self.exon_starts, self.exon_ends = None, [], []
        self.core, self.name2 = None, ''
        self.cds_start_stat, self.cds_end_stat = '', ''  # 4 choices: none, unk, incmpl, cmpl
        self.exon_frames = []
        self.read_count = 1
        self.read_density = 0

    @property
    def mRNA_length(self):
        return sum((np.array(self.exon_ends) - np.array(self.exon_starts) + 1))


class RefSeqTools(object):
    """
    We could download RefSeq Gene file from http://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/
    """
    def __init__(self):
        self.ref_genes = []

    def read_refgene_file(self, file_path):
        """
        read a refgene file which download from ENCODE
        """
        with open(file_path) as fi:
            for line in fi:
                content = line.strip().split()
                refgene = RefGene()
                refgene.bin = int(content[0])
                refgene.name, refgene.chrom, refgene.strand = content[1], content[2].lower(), content[3]
                refgene.tx_start, refgene.tx_end = int(content[4]), int(content[5])
                refgene.cds_start, refgene.cds_end = int(content[6]), int(content[7])
                refgene.exon_count = int(content[8])
                refgene.exon_starts = [int(e) for e in content[9][:-1].split(',')]
                refgene.exon_ends = [int(e) for e in content[10][:-1].split(',')]
                refgene.core, refgene.name2 = int(content[11]), content[12]
                refgene.cds_start_stat, refgene.cds_end_stat = content[13], content[14]
                # refgene.exon_frames = [int(e) for e in content[15][:-1].split(',')]
                self.ref_genes.append(refgene)

    def __get_reads_pos(self, reads_file):
        uniq_chrs = list(set([ref_gene.chrom for ref_gene in self.ref_genes]))
        read_starts = {key.lower(): [] for key in uniq_chrs}
        print 'reading %s reads ...' % reads_file
        file = open(reads_file)
        first_line_cnt = file.readline().strip().split()

        with open(reads_file) as fi:
            for line in fi:
                cnt = line.strip().split()
                chrm, start, end = cnt[0].lower(), int(cnt[1]), int(cnt[2])
                if chrm in uniq_chrs:
                    read_starts[chrm].append(start)

        for chrm in uniq_chrs:
            read_starts[chrm] = np.sort(np.array(read_starts[chrm]))
            print '%s: %d' % (chrm, read_starts[chrm].size)
        return read_starts

    def map_reads_2genes(self, reads_file):
        """
        map reads to ref_gene exons.
        """
        start1 = time()
        read_starts = self.__get_reads_pos(reads_file)
        start2 = time()
        times = 0
        for ref_gene in self.ref_genes:
            times += 1
            if times % 500 == 0:
                print 'calculated %d genes read count ...' % times
            if len(read_starts[ref_gene.chrom]) == 0:
                continue
            starts = read_starts[ref_gene.chrom]
            for es, ed in zip(ref_gene.exon_starts, ref_gene.exon_ends):
                # rd = starts[(starts > es) & (starts < ed)].size
                rd = cal_read_count(es, ed, starts)
                ref_gene.read_count += rd

        print 'start calculate rpkm ...'
        mapped_read_count = self.mapped_read_count
        for ref_gene in self.ref_genes:
            # calculate RPKM
            ref_gene.read_density = \
                ref_gene.read_count * 1000 * 1000 * 1000. / (ref_gene.mRNA_length * mapped_read_count)
        print 'got reads time: %f' % (time() - start1)
        print 'map reads time: %f' % (time() - start2)

    @property
    def mapped_read_count(self):
        return sum([gene.read_count for gene in self.ref_genes])


def cal_read_count(exon_start, exon_end, rnaseq_starts_position):
    si = bisect_left(rnaseq_starts_position, exon_start)
    ei = bisect_right(rnaseq_starts_position, exon_end)
    try:
        if exon_end == rnaseq_starts_position[ei]:
            return ei - si + 1
        else:
            return ei - si
    except IndexError:
        return ei - si


def output_refgenes_rpkm(refgene_file, name4save, rna_seq_files):
    """
    output gene expression value into a matrix txt file
    """
    file2save = open(name4save + '.txt', 'w')
    header = '\t'.join([''] + rna_seq_files) + '\n'
    file2save.write(header)
    result = []
    for rsf in rna_seq_files:
        ref_tool = RefSeqTools()
        ref_tool.read_refgene_file(refgene_file)
        ref_tool.map_reads_2genes(rsf)
        result.append(ref_tool)
    ref_genes = result[0].ref_genes
    for i, rg in enumerate(ref_genes):
        line = rg.name2 + '\t'
        line += '\t'.join([str(gene.read_density) for gene in [tool.ref_genes[i] for tool in result]])
        line += '\n'
        file2save.write(line)
    file2save.close()


def output_refgenes_read_count(refgene_file, name4save, rna_seq_files):
    """
    output gene expression value into a matrix txt file
    """
    file2save = open(name4save + '.txt', 'w')
    header = '\t'.join([''] + rna_seq_files) + '\n'
    file2save.write(header)
    result = []
    for rsf in rna_seq_files:
        ref_tool = RefSeqTools()
        ref_tool.read_refgene_file(refgene_file)
        ref_tool.map_reads_2genes(rsf)
        result.append(ref_tool)
    ref_genes = result[0].ref_genes
    for i, rg in enumerate(ref_genes):
        # line = rg.name2 + '\t'  # name2是gene symbol并非唯一的
        line = rg.name + '\t'  # name是transcript id是唯一的
        line += '\t'.join([str(gene.read_count) for gene in [tool.ref_genes[i] for tool in result]])
        line += '\n'
        file2save.write(line)
    file2save.close()


def run1(fi, refgene_file):
    if not os.path.exists(fi[:-4] + '_rpkm.txt'):
        output_refgenes_rpkm(refgene_file, fi[:-4] + '_rpkm', [fi])


def run2(fi, refgene_file):
    if not os.path.exists(fi[:-4] + '_read_count.txt'):
        output_refgenes_read_count(refgene_file, fi[:-4] + '_read_count', [fi])


if __name__ == '__main__':
    pass
