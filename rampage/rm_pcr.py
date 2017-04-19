#!/usr/bin/env python3

'''
Usage: rm_pcr.py [options] <rampage>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output prefix. [default: rampage]
'''

import os
import os.path
from seqlib.ngs import check_bam
from seqlib.seq import dna_to_rna

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def rm_pcr(options):
    '''
    Remove PCR duplicates
    '''
    # parse options
    rampage = options['<rampage>']
    thread = int(options['--thread'])
    prefix = options['--output']
    # remove PCR duplicates
    if thread > 1:  # more than 1 thread, parse bam file based on chromosome
        from multiprocessing import Pool
        p = Pool(thread)
        result = []
        for chrom in check_bam(rampage).references:
            if chrom.startswith('chr'):
                result.append(p.apply_async(remove_pcr, args=(rampage,
                                                              chrom,)))
        p.close()
        p.join()
        for n, r in enumerate(result):
            if n == 0:
                init = True
            else:
                init = False
            collapsed_pairs = r.get()
            write_signal(collapsed_pairs, prefix, init)
    else:
        collapsed_pairs = remove_pcr(rampage)
        write_signal(collapsed_pairs, prefix, True)


def check_file(fname, init):
    if os.path.isfile(fname) and init:
        os.remove(fname)
    return open(fname, 'a')


def write_signal(pairs, prefix, init):
    bed5p = check_file(prefix + '_plus_5end.bed', init)
    bed5m = check_file(prefix + '_minus_5end.bed', init)
    bed3p = check_file(prefix + '_plus_3read.bed', init)
    bed3m = check_file(prefix + '_minus_3read.bed', init)
    for pair in pairs:
        pair_info = pair.split()
        r1_chrom, r1_start, r1_end, r1_strand = pair_info[:4]
        r2_chrom, r2_start, r2_end, r2_strand = pair_info[5:9]
        if r1_strand == '+':
            bed5p.write('\t'.join([r1_chrom, r1_start, str(int(r1_start) + 1),
                                   '5end\t0', r1_strand]) + '\n')
            bed3p.write('\t'.join([r2_chrom, str(int(r2_end) - 1), r2_end,
                                   '3read\t0', r2_strand]) + '\n')
        else:
            bed5m.write('\t'.join([r1_chrom, str(int(r1_end) - 1), r1_end,
                                   '5end\t0', r1_strand]) + '\n')
            bed3m.write('\t'.join([r2_chrom, r2_start, str(int(r2_start) + 1),
                                   '3read\t0', r2_strand]) + '\n')


def remove_pcr(bam_f, chrom=None):
    bam = check_bam(bam_f)
    read1 = fetch_read1(bam, chrom)
    return fetch_read2(bam, read1, chrom)


def fetch_read1(bam, chrom):
    read1_lst = {}
    for read in bam.fetch(chrom, multiple_iterators=True):
        if read.is_read2 or read.is_secondary or read.mate_is_unmapped:
            continue
        chrom = read.reference_name
        mate_chrom = read.next_reference_name
        if chrom != mate_chrom:
            continue
        if not read.is_reverse and read.mate_is_reverse:
            strand = '+'
            mate_strand = '-'
        elif read.is_reverse and not read.mate_is_reverse:
            strand = '-'
            mate_strand = '+'
        else:
            continue
        mate_pos = str(read.next_reference_start)
        name = read.query_name
        start = str(read.reference_start)
        end = str(read.reference_end)
        if strand == '+':
            barcode = dna_to_rna(read.query_sequence[:6])
        else:
            barcode = dna_to_rna(read.query_sequence[-6:],
                                 strand=strand)
        read_id = '\t'.join([name, mate_chrom, mate_pos, mate_strand])
        read1_lst[read_id] = [name, chrom, start, end, strand, barcode]
    return read1_lst


def fetch_read2(bam, read1, chrom):
    collapsed_pairs = set()
    for read in bam.fetch(chrom, multiple_iterators=True):
        if read.is_read1 or read.is_secondary or read.mate_is_unmapped:
            continue
        name = read.query_name
        chrom = read.reference_name
        start = str(read.reference_start)
        strand = '+' if not read.is_reverse else '-'
        read_id = '\t'.join([name, chrom, start, strand])
        if read_id not in read1:
            continue
        end = str(read.reference_end)
        if strand == '+':
            barcode = dna_to_rna(read.query_sequence[:15])
        else:
            barcode = dna_to_rna(read.query_sequence[-15:],
                                 strand=strand)
        collapsed_pairs.add('\t'.join(read1[read_id][1:] + [chrom, start, end,
                                                            strand, barcode]))
    return collapsed_pairs


if __name__ == '__main__':
    from docopt import docopt
    rm_pcr(docopt(__doc__, version=__version__))
