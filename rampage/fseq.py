#!/usr/bin/env python3

'''
Usage: fseq.py [options] <rampagepeak>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -l LENGTH                      Feature length for F-seq. [default: 30]
    --wig                          Create Wig files.
'''

import sys
import os.path
import tempfile
import glob
from seqlib.path import which, check_dir
from seqlib.helper import run_command

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
    # run F-seq
    flist = ['rampage_plus_5end.bed', 'rampage_plus_3read.bed',
             'rampage_minus_5end.bed', 'rampage_minus_3read.bed']
    for f in flist:
        run_fseq(folder, f, flength, wig_flag)


def run_fseq(folder, f, flength, wig_flag):
    prefix = os.path.splitext(f)[0]
    # create bed files
    temp_dir = tempfile.mkdtemp()
    command = 'fseq -f 0 -l %s -of bed -o %s %s' % (flength, temp_dir,
                                                    os.path.join(folder, f))
    run_command(command, 'Error in F-seq!')
    cat_files(temp_dir, folder, prefix + '_fseq.bed')
    # create wig files
    if wig_flag:
        temp_dir = tempfile.mkdtemp()
        command = 'fseq -f 0 -l %s -o %s %s' % (flength, temp_dir,
                                                os.path.join(folder, f))
        run_command(command, 'Error in F-seq!')
        cat_files(temp_dir, folder, prefix + '_fseq.wig', is_wig=True)


def cat_files(temp_dir, folder, fname, is_wig=False):
    outf = os.path.join(folder, fname)
    with open(outf, 'w') as out:
        if is_wig:
            out.write('track type=wiggle_0 name=%s description=%s\n' % (fname,
                                                                        fname))
        for fname in glob.iglob(os.path.join(temp_dir, 'chr*')):
            with open(fname) as f:
                out.write(f.read())


if __name__ == '__main__':
    from docopt import docopt
    fseq(docopt(__doc__, version=__version__))
