#!/usr/bin/env python3

'''
Usage: assign_peak.py [options] (-r REF | -g GTF | --db=DB) <rampagepeak>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -r REF --ref=REF               Assembled gene annotation GenePred file.
    -g GTF --gtf=GTF               Assembled gene annotation GTF file.
    --db=DB                        Assembled gene annotation database.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    --promoter=PROMOTER            Promoter region. [default: 1000]
'''


import os.path
from multiprocessing import Pool
import numpy as np
import pysam
from seqlib.path import check_dir
from seqlib.ngs import check_bed
from interval import Interval

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def assign_peak(options):
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
    peak = check_bed(os.path.join(folder, 'rampage_peaks.txt'),
                     return_handle=False)
    prom = int(options['--promoter'])
    # align and filter candidate peak
    p = Pool(int(options['--thread']))
    results = []
    for gene_info, gpromoter in parse_gene(db, ref_flag, prom):
        results.append(p.apply_async(assign_peak_to_gene,
                                     args=(rampage, peak, gene_info,
                                           gpromoter, prom)))
    p.close()
    p.join()
    # output results
    with open(os.path.join(folder, 'rampage_assigned_peaks.txt'), 'w') as outf:
        for r in results:
            gene_info, peak = r.get()
            if gene_info:
                for p in peak:
                    outf.write('%s\t%s\n' % (p, gene_info))


def parse_gene(db, ref_flag, prom):
    if ref_flag:  # for GenePred
        gname = ''
        gchr, gstrand = '', ''
        ginterval = []  # genomic regions
        gpromoter = []  # tss regions
        with open(db, 'r') as f:
            for line in f:
                gene_id, chrom, strand, start, end = line.split()[:5]
                if not chrom.startswith('chr'):
                    continue
                start = int(start) + 1
                end = int(end)
                if gname == '':  # first entry
                    gname = gene_id
                    gchr, gstrand = chrom, strand
                # not same gene
                elif gname != gene_id or chrom != gchr or strand != gstrand:
                    gpromoter = Interval(gpromoter)  # combine tss regions
                    for itl in Interval(ginterval).interval:
                        gstart, gend = itl
                        gene_info = '%s\t%s\t%d\t%d\t%s' % (gname, gchr,
                                                            gstart, gend,
                                                            gstrand)
                        yield gene_info, gpromoter
                    # update gene info
                    gname = gene_id
                    gchr, gstrand = chrom, strand
                    ginterval = []
                    gpromoter = []
                # add genomic interval
                ginterval.append([start, end])
                # define gene promoter regions
                if strand == '+':
                    gpromoter.append([start - prom, start + prom, start])
                else:
                    gpromoter.append([end - prom, end + prom, end])
            else:  # last entry
                gpromoter = Interval(gpromoter)  # combine tss regions
                for itl in Interval(ginterval).interval:
                    gstart, gend = itl
                    gene_info = '%s\t%s\t%d\t%d\t%s' % (gname, gchr,
                                                        gstart, gend,
                                                        gstrand)
                    yield gene_info, gpromoter
    else:  # for GTF
        for gene in db.features_of_type('gene'):
            if not gene.seqid.startswith('chr'):
                continue
            gene_info = '%s\t%s\t%d\t%d\t%s' % (gene.id, gene.seqid,
                                                gene.start, gene.end,
                                                gene.strand)
            gpromoter = []  # tss regions
            for t in db.children(gene.id, featuretype='transcript'):
                if gene.strand == '+':
                    gpromoter.append([t.start - prom, t.start + prom,
                                      t.start])
                else:
                    gpromoter.append([t.end - prom, t.end + prom,
                                      t.end])
            gpromoter = Interval(gpromoter)  # combine tss regions
            yield gene_info, gpromoter


def assign_peak_to_gene(rampage, rampage_peak, gene_info, gp, prom):
    gene_id, gene_chrom, gene_start, gene_end, gene_strand = gene_info.split()
    gene_start = int(gene_start)
    gene_end = int(gene_end)
    peak_loc = report_peak_loc(rampage, gene_chrom, gene_start, gene_end,
                               gene_strand, prom)
    if not peak_loc:
        return None, None
    assigned_peaks = fetch_peak(rampage_peak, peak_loc, gene_chrom,
                                gene_strand, gp)
    return gene_info, assigned_peaks


def report_peak_loc(rampage, g_chrom, g_start, g_end, g_strand, prom):
    peak_loc = []
    rampagef = pysam.TabixFile(rampage)
    if g_chrom not in rampagef.contigs:
        return peak_loc
    for l in rampagef.fetch(g_chrom, g_start, g_end):
        start, end, _, _, strand = l.split()[1:6]
        if g_strand == '+':
            if strand == '-':  # not same strand
                continue
            end5 = int(start)
            read3 = int(end)
        else:
            if strand == '+':  # not same strand
                continue
            end5 = int(end)
            read3 = int(start)
        # ensure read3 within gene
        if read3 < g_start or read3 > g_end:
            continue
        # ensure end5 not far away from gene
        if end5 < g_start - prom or end5 > g_end + prom:
            continue
        peak_loc.append([end5 - 10, end5 + 10])
    return peak_loc


def fetch_peak(rampage_peak, peak_loc, chrom, strand, gp):
    peaks = set()
    peakf = pysam.TabixFile(rampage_peak)
    for loc in Interval(peak_loc):
        start, end = loc[:2]
        for p in peakf.fetch(chrom, start, end):
            s = p.split()[5]
            if s != strand:  # not same strand
                continue
            pos = int(p.split()[6])
            tss_info = fetch_tss(gp.interval, pos)
            peaks.add(p.rstrip() + '\t%s' % tss_info)
    return peaks


def fetch_tss(gpromoter, pos):
    promoter_list = Interval.mapto([pos, pos], gpromoter)
    if not promoter_list:
        return 'None\tNone'
    else:
        promoter_list = promoter_list[0][2:]
    nearest_p = ''
    nearest_d = np.inf
    for ploc in promoter_list:
        distance = abs(pos - ploc)
        if distance < nearest_d:
            nearest_p = ploc
            nearest_d = distance
    return '%d\t%d' % (nearest_p, nearest_d)


if __name__ == '__main__':
    from docopt import docopt
    assign_peak(docopt(__doc__, version=__version__))
