#!/usr/bin/env python3

'''
Usage: call_peak.py [options] (-r REF | -g GTF | --db=DB) <rampagepeak>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -r REF --ref=REF               Assembled gene annotation GenePred file.
    -g GTF --gtf=GTF               Assembled gene annotation GTF file.
    --db=DB                        Assembled gene annotation database.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    --promoter=PROMOTER            Promoter region. [default: 1000]
'''

import math
import os.path
from multiprocessing import Pool
import numpy as np
import scipy.stats
import pysam
from seqlib.path import check_dir
from seqlib.ngs import check_bed
from interval import Interval

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def call_peak(options):
    '''
    Call rampage peaks
    '''
    # parse options
    if options['--ref']:
        db = options['--ref']
        ref_flag = True
    elif options['--db']:
        import gffutils
        db = gffutils.FeatureDB(options['--db'])
        ref_flag = False
    else:
        import gffutils
        gtf_f = options['--gtf']
        prefix = os.path.splitext(os.path.basename(gtf_f))[0]
        db = gffutils.create_db(gtf_f, prefix + '.db',
                                force=True, disable_infer_transcripts=True)
        ref_flag = False
    folder = check_dir(options['<rampagepeak>'])
    rampage = check_bed(os.path.join(folder, 'rampage_link.bed'),
                        return_handle=False)
    peak5 = {'+': check_bed(os.path.join(folder, 'rampage_plus_5end_fseq.bed'),
                            return_handle=False),
             '-': check_bed(os.path.join(folder,
                                         'rampage_minus_5end_fseq.bed'),
                            return_handle=False)}
    peak3 = {'+': check_bed(os.path.join(folder,
                                         'rampage_plus_3read_fseq.bed'),
                            return_handle=False),
             '-': check_bed(os.path.join(folder,
                                         'rampage_minus_3read_fseq.bed'),
                            return_handle=False)}
    pregion = int(options['--promoter'])
    # align and filter candidate peak
    p = Pool(int(options['--thread']))
    results = []
    for gene_info, gpromoter, gstrand in parse_gene(db, ref_flag, pregion):
        p5, p3 = peak5[gstrand], peak3[gstrand]
        results.append(p.apply_async(call_peak_for_gene,
                                     args=(rampage, p5, p3, gene_info,
                                           gpromoter, pregion,)))
    p.close()
    p.join()
    peaks, pvalue = [], []
    for r in results:
        gene_info, peak = r.get()
        if gene_info:
            for p in peak:
                peaks.append(gene_info + '\t' + p[0])
                pvalue.append(p[1])
    # calculate q-value
    rank = {p: r for r, p in enumerate(np.array(pvalue).argsort())}
    total_num = len(pvalue)
    qvalue = [p * total_num / (rank[n] + 1) for n, p in enumerate(pvalue)]
    # output results
    with open(os.path.join(folder, 'rampage_peak.txt'), 'w') as outf:
        for p, q in zip(peaks, qvalue):
            outf.write('%s\t%f\n' % (p, q))


def parse_gene(db, ref_flag, pregion):
    if ref_flag:  # for GenePred
        gname = ''
        gstart, gend, gchr, gstrand = 0, 0, '', ''
        gpromoter = []
        with open(db, 'r') as f:
            for line in f:
                gene_id, chrom, strand, start, end = line.split()[:5]
                if not chrom.startswith('chr'):
                    continue
                start = int(start) + 1
                end = int(end)
                if gname == gene_id or gname == '':  # gene not changed
                    if gname == '':  # first entry
                        gname = gene_id
                        gstart, gend, gchr, gstrand = start, end, chrom, strand
                    else:  # not first entry
                        gstart = start if start < gstart else gstart
                        gend = end if end > gend else gend
                else:  # gene changed
                    gene_info = '%s\t%s\t%d\t%d\t%s' % (gname, gchr, gstart,
                                                        gend, gstrand)
                    gpromoter = Interval(gpromoter)
                    yield gene_info, gpromoter, gstrand
                    gname = gene_id
                    gstart, gend, gchr, gstrand = start, end, chrom, strand
                    gpromoter = []
                # define gene promoter regions
                FPKM, TPM = line.split()[-2:]
                if strand == '+':
                    pinfo = '%d|%s|%s' % (start, FPKM, TPM)
                    gpromoter.append([start - pregion, start + pregion, pinfo])
                else:
                    pinfo = '%d|%s|%s' % (end, FPKM, TPM)
                    gpromoter.append([end - pregion, end + pregion, pinfo])
            else:  # last entry
                gene_info = '%s\t%s\t%d\t%d\t%s' % (gname, gchr, gstart,
                                                    gend, gstrand)
                gpromoter = Interval(gpromoter)
                yield gene_info, gpromoter, gstrand
    else:  # for GTF
        for gene in db.features_of_type('gene'):
            if not gene.seqid.startswith('chr'):
                continue
            gene_info = '%s\t%s\t%d\t%d\t%s' % (gene.id, gene.seqid,
                                                gene.start, gene.end,
                                                gene.strand)
            gpromoter = []
            for t in db.children(gene.id, featuretype='transcript'):
                FPKM = t.attributes['FPKM']
                TPM = t.attributes['TPM']
                if gene.strand == '+':
                    pinfo = '%d|%s|%s' % (t.start, FPKM, TPM)
                    gpromoter.append([t.start - pregion, t.start + pregion,
                                      pinfo])
                else:
                    pinfo = '%d|%s|%s' % (t.end, FPKM, TPM)
                    gpromoter.append([t.end - pregion, t.end + pregion,
                                      pinfo])
            gpromoter = Interval(gpromoter)
            yield gene_info, gpromoter, gstrand


def call_peak_for_gene(rampage, peak5, peak3, gene_info, gpromoter, pregion):
    gene_id, gene_chrom, gene_start, gene_end, gene_strand = gene_info.split()
    gene_start = int(gene_start)
    gene_end = int(gene_end)
    peak_loc = assign_peak(rampage, gene_chrom, gene_start, gene_end,
                           gene_strand, pregion)
    peaks = fetch_peak(peak5, peak_loc, gene_chrom)
    e_mean, var = cal_expression(peak3, gene_chrom, gene_start, gene_end)
    if e_mean == 0:
        return None, None
    if var == 0:
        return None, None
    filtered_peaks = filter_peak(peaks, e_mean, var, gpromoter)
    return gene_info, filtered_peaks


def assign_peak(rampage, g_chrom, g_start, g_end, g_strand, pregion):
    peak_loc = []
    rampagef = pysam.TabixFile(rampage)
    for l in rampagef.fetch(g_chrom, g_start, g_end):
        start, end, _, _, strand = l.split()[1:6]
        if g_strand == '+':
            if strand == '-':
                continue
            end5 = int(start)
            read3 = int(end)
        else:
            if strand == '+':
                continue
            end5 = int(end)
            read3 = int(start)
        # ensure read3 within gene
        if read3 < g_start or read3 > g_end:
            continue
        # ensure end5 not far away from gene
        if end5 < g_start - pregion or end5 > g_end + pregion:
            continue
        peak_loc.append([end5 - 10, end5 + 10])
    return peak_loc


def fetch_peak(peak5, peak_loc, chrom):
    peaks = set()
    peak5f = pysam.TabixFile(peak5)
    for loc in Interval(peak_loc):
        start, end = loc[:2]
        for p in peak5f.fetch(chrom, start, end):
            peaks.add(p)
    return peaks


def cal_expression(peak3, g_chrom, g_start, g_end):
    expression = []
    peak3f = pysam.TabixFile(peak3)
    for e in peak3f.fetch(g_chrom, g_start, g_end):
        expression.append(float(e.split()[4]))
    expression = np.array(expression)
    if expression.size == 0:
        return 0, 0
    e_mean = np.mean(expression)
    var = math.sqrt(np.var(expression) / expression.size)
    return e_mean, var


def filter_peak(peaks, e_mean, var, gpromoter):
    filtered_peaks = []
    for p in peaks:
        chrom, start, end, _, value = p.split()
        value = float(value)
        z_score, p_value = wald_test(value, e_mean, var)
        pos = int((int(start) + int(end)) / 2)
        einfo = fetch_expression(gpromoter.interval, pos)
        filtered_peaks.append(['\t'.join([chrom, start, end, einfo, str(value),
                                          str(e_mean), str(z_score),
                                          str(p_value)]),
                               p_value])
    return filtered_peaks


def wald_test(value, mean, var):
    '''
    Test whether a value belongs to samples
    '''
    z = (value - mean) / var
    p = scipy.stats.norm.sf(z)
    return z, p


def fetch_expression(gpromoter, pos):
    promoter_FPKM = {}
    promoter_TPM = {}
    promoter_list = Interval.mapto([pos, pos], gpromoter)
    if not promoter_list:
        return 'None'
    else:
        promoter_list = promoter_list[0][2:]
    for p in promoter_list:
        ploc, FPKM, TPM = p.split('|')
        ploc = int(ploc)
        if ploc in promoter_FPKM:
            if FPKM > promoter_FPKM[ploc]:
                promoter_FPKM[ploc] = FPKM
                promoter_TPM[ploc] = TPM
        else:
            promoter_FPKM[ploc] = FPKM
            promoter_TPM[ploc] = TPM
    nearest_p = ''
    nearest_d = np.inf
    for ploc in promoter_FPKM:
        distance = abs(pos - ploc)
        if distance < nearest_d:
            nearest_p = '%d|%s|%s' % (ploc, promoter_FPKM[ploc],
                                      promoter_TPM[ploc])
            nearest_d = distance
    return nearest_p


if __name__ == '__main__':
    from docopt import docopt
    call_peak(docopt(__doc__, version=__version__))
