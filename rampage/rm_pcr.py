#!/usr/bin/env python3

'''
Usage: rm_pcr.py [options] <rampage>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output directory. [default: rampage_peak]
    --min=MIN                      Minimum read counts. [default: 1]
'''

import os
import os.path
from collections import defaultdict
from seqlib.path import create_dir
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
    folder = create_dir(options['--output'])
    min_read = int(options['--min'])
    total_counts = 0
    # remove PCR duplicates
    if thread > 1:  # more than 1 thread, parse bam file based on chromosome
        from multiprocessing import Pool
        p = Pool(thread)
        result = []
        for chrom in check_bam(rampage).references:  # for each chromosome
            if chrom.startswith('chr'):  # only parse normal chromosome
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
            total_counts += write_signal(collapsed_pairs, folder, min_read,
                                         init)
    else:
        collapsed_pairs = remove_pcr(rampage)
        total_counts += write_signal(collapsed_pairs, folder, min_read, True)
    print('Total counts: %d' % total_counts)


def check_file(fname, init):
    if os.path.isfile(fname) and init:
        os.remove(fname)
    return open(fname, 'a')


def write_signal(pairs, folder, min_read, init):
    bed5p = check_file(folder + '/rampage_plus_5end.bed', init)
    bed5m = check_file(folder + '/rampage_minus_5end.bed', init)
    bed3p = check_file(folder + '/rampage_plus_3read.bed', init)
    bed3m = check_file(folder + '/rampage_minus_3read.bed', init)
    bed = check_file(folder + '/rampage_link.bed', init)
    uniq_pairs = defaultdict(int)
    for pair in pairs:
        info = pair.rsplit('\t', 1)[0]
        uniq_pairs[info] += 1
    total_counts = 0
    for pair in uniq_pairs:
        read_count = uniq_pairs[pair]
        if read_count < min_read:  # not enough reads
            continue
        total_counts += read_count
        pair_info = pair.split()
        r1_chrom, r1_start, r1_end, strand = pair_info[:4]
        r2_chrom, r2_start, r2_end = pair_info[4:]
        if strand == '+':
            start = int(r1_start)
            end = int(r2_end)
            bed5p.write('%s\t%d\t%d\t5end\t0\t%s\n' % (r1_chrom, start,
                                                       start + 1, strand) *
                        read_count)
            bed3p.write('%s\t%d\t%d\t3read\t0\t%s\n' % (r2_chrom, end - 1,
                                                        end, strand) *
                        read_count)
            offset = '0,' + str(end - start - 1)
            bed.write('\t'.join([r1_chrom, r1_start, r2_end, 'link\t0',
                                 strand, r1_start, r1_start, '0,0,0',
                                 '2', '1,1', offset]) + '\n')
        else:
            start = int(r2_start)
            end = int(r1_end)
            bed5m.write('%s\t%d\t%d\t5end\t0\t%s\n' % (r1_chrom, end - 1,
                                                       end, strand) *
                        read_count)
            bed3m.write('%s\t%d\t%d\t3read\t0\t%s\n' % (r2_chrom, start,
                                                        start - 1, strand) *
                        read_count)
            offset = '0,' + str(end - start - 1)
            bed.write('\t'.join([r1_chrom, r2_start, r1_end, 'link\t0',
                                 strand, r2_start, r2_start, '0,0,0',
                                 '2', '1,1', offset]) + '\n')
    return total_counts


def remove_pcr(bam_f, chrom=None):
    bam = check_bam(bam_f)
    read1 = fetch_read1(bam, chrom)  # fetch read1
    return fetch_read2(bam, read1, chrom)  # fetch read2


def fetch_read1(bam, chrom):
    read1_lst = {}
    for read in bam.fetch(chrom):
        # not read2 or secondary alignment or read2 unmapped
        if read.is_read2 or read.is_secondary or read.mate_is_unmapped:
            continue
        chrom = read.reference_name
        mate_chrom = read.next_reference_name
        if chrom != mate_chrom:  # not same chromosome
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
        read_id = '\t'.join([name, mate_chrom, mate_pos, mate_strand])
        read1_lst[read_id] = [name, chrom, start, end, strand]
    return read1_lst


def fetch_read2(bam, read1, chrom):
    collapsed_pairs = set()
    # not read1 or secondary alignment or read1 unmapped
    for read in bam.fetch(chrom):
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
                                                            barcode]))
    return collapsed_pairs


if __name__ == '__main__':
    from docopt import docopt
    rm_pcr(docopt(__doc__, version=__version__))
