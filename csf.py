#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Usage: csf.py [options] -g GENOME -c CHROM -s SITE -r STRAND (-1 FA1 -2 FA2 | \
-R FA) <out_dir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -g GENOME --genome=GENOME      Genome FASTA file.
    -c CHROM --chrom=CHROM         Chromosome of sgRNA.
    -s SITE --site=SITE            Site of sgRNA.
    -r STRAND --strand=STRAND      Strand of sgRNA (1: +, 0: -).
    -1 FA1                         Comma-separated list of mate 1 files.
    -2 FA2                         Comma-separated list of mate 2 files.
    -R FA                          Comma-separated list of fastq files.
    -p THREAD --thread=THREAD      Running threads. [default: 10]
    --stranded                     Stranded sequencing tag. # TODO
    --read-length=RLEN             Length of reads. [default: 100]
    --region-length=ALEN           Length of sgRNA region. [default: 50]
    --check-region-length=CLEN     Length of check region. [default: 20]
    --skip-alignment               Skip alignment step.
'''

import sys
import os
import os.path
import tempfile
import gzip
import pysam
from collections import defaultdict
from docopt import docopt
from seqlib.path import create_dir, check_dir, which
from seqlib.ngs import check_fasta

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def main():
    # parse options
    options = docopt(__doc__, version=__version__)
    fa = check_fasta(options['--genome'])
    chrom = options['--chrom']
    site = int(options['--site'])
    strand = '+' if options['--strand'] == '1' else '-'
    rlen = int(options['--read-length'])
    alen = int(options['--region-length'])
    clen = int(options['--check-region-length'])
    thread = options['--thread']
    skip_flag = options['--skip-alignment']
    # check output directory
    if not skip_flag:  # not skip alignment
        out_dir = create_dir(options['<out_dir>'])
    else:  # skip alignment
        out_dir = check_dir(options['<out_dir>'])
    # build index for sgRNA
    index_path, offset = build_index(fa, chrom, site, strand, rlen, thread,
                                     out_dir)
    if not skip_flag:  # not skip alignment
        # deal with reads file
        reads = tempfile.NamedTemporaryFile(mode='w+')
        if options['-R']:
            fq_lst = options['-R'].split(',')
            convert_read(reads, single=fq_lst)
        else:
            fq1_lst = options['-1'].split(',')
            fq2_lst = options['-2'].split(',')
            convert_read(reads, fq1=fq1_lst, fq2=fq2_lst)
        reads.seek(0)
        read_path = reads.name
        # mapped reads with bowtie2
        bam = bowtie2_align(index_path, read_path, thread, out_dir)
        # remove tempfile
        reads.close()
    else:
        bam = os.path.join(out_dir, 'cs.bam')
    # fetch cleavage site reads
    fetch_reads(index_path, offset, strand, alen, clen, bam, out_dir)


def build_index(fa, chrom, site, strand, rlen, thread, out_dir):
    print('Build index...')
    if strand == '+':
        start = site - (rlen - 10)
        end = site + (rlen - 20)
        offset = rlen - 10
    else:
        start = site - (rlen - 20)
        end = site + (rlen - 10)
        offset = rlen - 20
    index_path = os.path.join(out_dir, 'sgRNA.fa')
    # fetch sgRNA region sequence
    with open(index_path, 'w') as out:
        out.write('>sgRNA_region\n')
        out.write(fa.fetch(chrom, start, end) + '\n')
    # build index
    if which('bowtie2-build'):
        command = 'bowtie2-build -q --threads %s %s %s'
        return_code = os.system(command % (thread, index_path,
                                           index_path)) >> 8
        if return_code:
            sys.exit('Error: cannot build index for sgRNA!')
    else:
        sys.exit('Error: no bowtie2-build installed!')
    return index_path, offset


def convert_read(reads, single=None, fq1=None, fq2=None):
    print('Convert read files...')
    if single:  # single-end
        write_read(reads, single)
    else:  # paired-end
        write_read(reads, fq1, tag='first')
        write_read(reads, fq2, tag='second')


def write_read(reads, lst, tag='none'):
    if tag == 'none':  # single-end
        first_tag = False
        tag_flag = False
    elif tag == 'first':  # first mate
        first_tag = True
        tag_flag = True
    else:  # second mate
        first_tag = False
        tag_flag = True
    for read_f in lst:
        if os.path.isfile(read_f):  # file exists
            if read_f.endswith('gz'):  # file in gzip format
                with gzip.open(read_f, 'rb') as f:
                    add_tag(reads, f, first_tag, tag_flag)
            else:  # file in text format
                with open(read_f, 'r') as f:
                    add_tag(reads, f, first_tag, tag_flag)
        else:  # file does not exist
            sys.exit('No fastq file: %s' % read_f)


def add_tag(reads, f, first_tag, tag_flag):
    if tag_flag:  # add tag
        for n, line in enumerate(f):
            if n % 4 == 0:  # read id line
                read_id = line.split()[0]
                if first_tag:  # first mate
                    reads.write(read_id + '/1\n')
                else:  # second mate
                    reads.write(read_id + '/2\n')
            else:  # not read id line
                reads.write(line)
    else:  # not add tag
        reads.write(f.read())


def bowtie2_align(index, read, thread, out_dir):
    print('Align reads...')
    if which('bowtie2'):
        bam = os.path.join(out_dir, 'cs.bam')
        sam = tempfile.NamedTemporaryFile('w+')
        command = 'bowtie2 --quiet --end-to-end -p %s -x %s -U %s -S %s'
        return_code = os.system(command % (thread, index, read, sam.name)) >> 8
        if return_code:
            sys.exit('Error in bowtie2 alignment!')
        sam.seek(0)
        with pysam.AlignmentFile(sam.name, 'r') as sam_f:
            with pysam.AlignmentFile(bam, 'wb', template=sam_f) as bam_f:
                for read in sam_f:
                    if not read.is_unmapped:
                        bam_f.write(read)
        sam.close()
        return bam
    else:
        sys.exit('Error: no bowtie2 installed!')


def fetch_reads(index, offset, strand, alen, clen, bam, out_dir):
    cs = os.path.join(out_dir, 'cs_region.txt')
    count = os.path.join(out_dir, 'cs_count.txt')
    with open(cs, 'w') as cs_out, open(count, 'w') as count_out:
        if strand == '+':
            cs_out.write(' ' * 35 + 'sgRNA|PAM\n')
        else:
            cs_out.write(' ' * 37 + 'PAM|sgRNA\n')
        index_fa = check_fasta(index).fetch('sgRNA_region')
        cs_out.write('Reference: ' + index_fa[offset - alen:offset + alen] +
                     '\n')
        count_out.write('Reference: ' + index_fa[offset - clen:offset + clen] +
                        '\n')
        cs_count = defaultdict(int)
        bam_f = pysam.AlignmentFile(bam, 'rb')
        for read in bam_f:
            seq = read.query_sequence
            start = read.reference_start
            pos = 0
            align = insert = ' ' * start
            for tag, tlen in read.cigartuples:
                if tag == 0:  # M
                    align += seq[pos:pos+tlen]
                    insert += ' ' * tlen
                    pos += tlen
                elif tag == 1:  # I
                    insert += seq[pos:pos+tlen]
                    pos += tlen
                else:  # D
                    align += '*' * tlen
            cs_out.write('Alignment: ' + align[offset - alen:offset + alen] +
                         '\n')
            cs_out.write('Insertion: ' + insert[offset - alen:offset + alen] +
                         '\n')
            cs_id = 'Reads    : ' + align[offset - clen:offset + clen]
            cs_id += '\t%d\n'
            cs_id += 'Indel    :' + insert[offset - clen:offset + clen] + '\n'
            cs_count[cs_id] += 1
        for cs_id in cs_count:
            count_out.write(cs_id % cs_count[cs_id])


if __name__ == '__main__':
    main()
