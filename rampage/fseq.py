#!/usr/bin/env python3

'''
Usage: fseq.py [options] <rampagepeak>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -l LENGTH                      Feature length for F-seq. [default: 30]
    --wig                          Create Wig files.
    -p PERCENT                     Retained percent of reads in resized peaks.
                                   [default: 0.95]
'''

import sys
import os.path
import tempfile
import glob
import operator
from collections import Counter
from seqlib.path import which, check_dir
from seqlib.helper import run_command
from seqlib.ngs import check_bed

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def fseq(options):
    '''
    Call peaks using F-seq
    '''
    # parse options
    if not which('fseq'):
        sys.exit('Error: No F-seq installed!')
    folder = check_dir(options['<rampagepeak>'])
    flength = options['-l']
    wig_flag = options['--wig']
    percent = float(options['-p'])
    # run F-seq
    flist = {'+': 'rampage_plus_5end.bed', '-': 'rampage_minus_5end.bed'}
    for strand in flist:
        run_fseq(folder, flist[strand], strand, flength, wig_flag, percent)


def run_fseq(folder, bed, strand, flength, wig_flag, percent):
    prefix = os.path.splitext(bed)[0]
    # create bed files
    temp_dir = tempfile.mkdtemp()
    bed_f = os.path.join(folder, bed)
    # command = 'fseq -f 0 -l %s -of bed -o %s %s' % (flength, temp_dir, bed_f)
    # run_command(command, 'Error in F-seq!')
    peak = prefix + '_fseq.bed'
    # cat_files(temp_dir, folder, peak)
    resized_peak = prefix + '_peak.bed'
    resize_peak(folder, peak, bed_f, resized_peak, strand, percent)
    # create wig files
    if wig_flag:
        temp_dir = tempfile.mkdtemp()
        command = 'fseq -f 0 -l %s -o %s %s' % (flength, temp_dir, bed_f)
        run_command(command, 'Error in F-seq!')
        wig = prefix + '_fseq.wig'
        cat_files(temp_dir, folder, wig, is_wig=True)


def cat_files(temp_dir, folder, fname, is_wig=False):
    outf = os.path.join(folder, fname)
    with open(outf, 'w') as out:
        if is_wig:
            out.write('track type=wiggle_0 name=%s description=%s\n' % (fname,
                                                                        fname))
        for fname in glob.iglob(os.path.join(temp_dir, 'chr*')):
            with open(fname) as f:
                out.write(f.read())


def resize_peak(folder, peak, bed, resized_peak, strand, percent):
    bed_f = check_bed(bed)
    peak_f = os.path.join(folder, peak)
    resized_peak_f = os.path.join(folder, resized_peak)
    with open(peak_f, 'r') as f, open(resized_peak_f, 'w') as out:
        for line in f:
            chrom, start, end = line.split()[:3]
            start = int(start)
            end = int(end) + 1
            total = 0
            sites = []
            for read in bed_f.fetch(chrom, start, end):
                sites.append(int(read.split()[1]))
                total += 1
            sites = Counter(sites)
            loc, height = sites.most_common()[0]
            peak_sites = 0
            for i in range(loc - 2, loc + 2):
                peak_sites += sites.get(i, 0)
            required_sites = int(total * percent)
            sub_sites = 0
            region = []
            for site, num in sorted(sites.items(), key=operator.itemgetter(1),
                                    reverse=True):
                sub_sites += num
                region.append(site)
                if sub_sites >= required_sites:
                    break
            region.sort()
            new_start, new_end = region[0], region[-1] + 1
            out_format = '%s\t%d\t%d\tpeak\t0\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n'
            out.write(out_format % (chrom, new_start, new_end, strand, loc,
                                    height, peak_sites, total, start, end))


if __name__ == '__main__':
    from docopt import docopt
    fseq(docopt(__doc__, version=__version__))
