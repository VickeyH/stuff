#!/usr/bin/env python

'''
Usage: count_site.py [options] -r REF -g GENOME -b BIGWIG -c COUNT -s SCORE

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -r REF --ref=REF               Gene annotation.
    -g GENOME --genome=GENOME      Genome FASTA file.
    -b BIGWIG --bigwig=BIGWIG      PhastCons BigWig file.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    --intron-length=LENGTH         Minimum intron length. [default: 10000]
    --min-distance=DISTANCE        Minimum distance to splicing sites.
                                   [default: 2000]
    -c COUNT
    -s SCORE
'''

import re
from multiprocessing import Pool
import pyBigWig
from seqlib.ngs import check_fasta
from seqlib.parse import Annotation
from seqlib.seq import dna_to_rna
from itertools import product
from collections import defaultdict
from collections import Counter

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def find_rs(options):
    # parse options
    ref = Annotation(options['--ref'])
    min_intron_length = int(options['--intron-length'])
    # prepare to parse intron regions
    rs_intron_set = set()
    p = Pool(int(options['--thread']))
    results = []
    for isoform in ref:
        chrom, gene, strand = isoform.chrom, isoform.gene, isoform.strand
        for start, end in zip(isoform.intron_starts, isoform.intron_ends):
            intron_length = end - start
            if intron_length < min_intron_length:  # short intron
                continue
            intron_info = '\t'.join([gene, chrom, str(start), str(end), strand,
                                     str(intron_length)])
            if intron_info in rs_intron_set:  # duplicated intron
                continue
            rs_intron_set.add(intron_info)
            results.append(p.apply_async(parse_intron, args=(options, chrom,
                                                             start, end,
                                                             strand,
                                                             intron_info)))
    p.close()
    p.join()
    # output results
    with open(options['-c'], 'w') as count_f, open(options['-s'], 'w') as score_f:
        for n, r in enumerate(results):
            rs_list, intron_length = r.get()
            if n == 0:
                total_rs_list = Counter(rs_list)
                total_intron_length = intron_length
            else:
                total_rs_list += Counter(rs_list)
                total_intron_length += intron_length
        for mer in total_rs_list:
            count_f.write('%s\t%f\n' % (mer, (len(total_rs_list[mer]) - n - 1) * 1.0 / total_intron_length))
            score_f.write(mer + '\t' + '\t'.join(str(x) for x in total_rs_list[mer] if x != 'x') + '\n')


def parse_intron(options, chrom, start, end, strand, intron_info):
    # fetch fasta
    fa = check_fasta(options['--genome'])
    intron_fa = dna_to_rna(fa.fetch(chrom, start, end), strand)
    # parse options
    phastcons_f = pyBigWig.open(options['--bigwig'])
    min_distance = int(options['--min-distance'])
    # start to parse rs sites
    rs_list = defaultdict(list)
    for i in product('ATCG', repeat=4):
        motif = ''.join(i)
        rs_list[motif].append('x')
        for m in re.finditer(motif, intron_fa):
            if strand == '+':
                pos = start + m.start() + 2
                left_dist, right_dist, dist_flag = cal_distance(pos, start,
                                                                end,
                                                                min_distance)
                if not dist_flag:  # not enough distance
                    continue
            else:
                pos = end - m.start() - 2
                left_dist, right_dist, dist_flag = cal_distance(pos, start,
                                                                end,
                                                                min_distance)
                if not dist_flag:  # not enough distance
                    continue
            phastcons = phastcons_f.stats(chrom, pos - 2, pos + 2)[0]
            if phastcons is None:  # no conservation score
                continue
            rs_list[motif].append(phastcons)
    intron_length = end - start - 2 * min_distance
    return(rs_list, intron_length)


def cal_distance(pos, start, end, min_distance):
    left_distance = pos - start
    if left_distance < min_distance:
        return (None, None, False)
    right_distance = end - pos
    if right_distance < min_distance:
        return (None, None, False)
    return (left_distance, right_distance, True)


if __name__ == '__main__':
    from docopt import docopt
    find_rs(docopt(__doc__, version=__version__))
