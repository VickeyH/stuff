#!/usr/bin/env python3

'''
Usage: call_peak.py [options] (-g GTF | --db DB) <rampagepeak>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -g GTF --gtf=GTF               Assembled gene annotation GTF file.
    --db DB                        Assembled gene annotation database.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    --cutoff=CUFOFF                P-value cutoff. [default: 0.05]
'''

import math
import os.path
import numpy as np
import scipy.stats
import gffutils
from seqlib.path import check_dir
from seqlib.ngs import check_bed
from interval import Interval

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def call_peak(options):
    '''
    Call rampage peaks
    TODO: 1. cutoff; 2. multiprocessing; 3. FDR
    '''
    # parse options
    if options['--gtf']:
        gtf_f = options['--gtf']
        prefix = os.path.splitext(os.path.basename(gtf_f))[0]
        db = gffutils.create_db(gtf_f, prefix + '.db',
                                force=True, disable_infer_transcripts=True)
    else:
        db = gffutils.FeatureDB(options['--db'])
    folder = check_dir(options['<rampagepeak>'])
    rampage = check_bed(os.path.join(folder, 'rampage_link.bed'))
    peak5 = {'+': check_bed(os.path.join(folder,
                                         'rampage_plus_5end_fseq.bed')),
             '-': check_bed(os.path.join(folder,
                                         'rampage_minus_5end_fseq.bed'))}
    peak3 = {'+': check_bed(os.path.join(folder,
                                         'rampage_plus_3read_fseq.bed')),
             '-': check_bed(os.path.join(folder,
                                         'rampage_minus_3read_fseq.bed'))}
    cutoff = float(options['--cutoff'])
    # align and filter candidate peak
    outf = open(os.path.join(folder, 'rampage_peak.txt'), 'w')
    for gene in db.features_of_type('gene'):
        peak_loc = assign_peak(rampage, gene)
        peaks = fetch_peak(peak5[gene.strand], peak_loc, gene)
        e_mean, var = cal_expression(peak3[gene.strand], gene)
        if e_mean == 0:
            continue
        if var == 0:
            continue
        filtered_peaks = filter_peak(peaks, e_mean, var, cutoff)
        gene_info = '%s\t%d\t%d\t%f' % (gene.seqid, gene.start, gene.end,
                                        e_mean)
        for p in filtered_peaks:  # TODO: support multiprocessing
            outf.write(p + '\t' + gene_info + '\n')


def assign_peak(rampage, gene):
    peak_loc = []
    for l in rampage.fetch(gene.seqid, gene.start, gene.end):
        start, end, _, _, strand = l.split()[1:6]
        start = int(start)
        end = int(end)
        if start < gene.start or gene.end < end:
            continue
        if gene.strand == '+':
            if strand == '-':
                continue
            end5 = start
        else:
            if strand == '+':
                continue
            end5 = end
        peak_loc.append([end5 - 10, end5 + 10])
    return peak_loc


def fetch_peak(peak5, peak_loc, gene):
    peaks = set()
    for loc in Interval(peak_loc):
        start, end = loc[:2]
        for p in peak5.fetch(gene.seqid, start, end):
            peaks.add(p)
    return peaks


def cal_expression(peak3, gene):
    expression = []
    for e in peak3.fetch(gene.seqid, gene.start, gene.end):
        expression.append(float(e.split()[4]))
    expression = np.array(expression)
    if expression.size == 0:
        return 0, 0
    e_mean = np.mean(expression)
    var = math.sqrt(np.var(expression) / expression.size)
    return e_mean, var


def filter_peak(peaks, e_mean, var, cutoff):
    filtered_peaks = []
    for p in peaks:
        chrom, start, end, _, value = p.split()
        value = float(value)
        z_score, p_value = wald_test(value, e_mean, var)
        filtered_peaks.append('\t'.join([chrom, start, end, str(value),
                                         str(z_score), str(p_value)]))
    return filtered_peaks


def wald_test(value, mean, var):
    '''
    Test whether a value belongs to samples
    '''
    z = (value - mean) / var
    p = scipy.stats.norm.sf(z)
    return z, p


if __name__ == '__main__':
    from docopt import docopt
    call_peak(docopt(__doc__, version=__version__))
