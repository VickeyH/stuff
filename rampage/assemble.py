#!/usr/bin/env python3

'''
Usage: assemble.py [options] -g GTF <bam>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -g GTF --gtf=GTF               Gene annotation GTF file.
    -p THREAD --thread=THREAD      Threads. [default: 5]
'''

import sys
import os
import os.path
from seqlib.path import which
from seqlib.helper import run_command

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def assemble(options):
    '''
    Assemble RNA-seq with StringTie
    TODO: convert GTF to GenePred
    TODO: check StringTie command
    '''
    # parse options
    if not which('stringtie'):
        sys.exit('Error: No StringTie installed!')
    bam = options['<bam>']
    gtf = options['--gtf']
    thread = options('--thread')
    run_stringtie(bam, gtf, thread)


def run_stringtie(bam, gtf, thread):
    if not os.path.isfile(bam):
        sys.exit('No BAM file: %s!' % bam)
    if not os.path.isfile(gtf):
        sys.exit('No GTF file: %s' % gtf)
    dirname, fname = os.path.split(bam)
    if not dirname:
        dirname = os.getcwd()
    prefix = os.path.splitext(fname)[0]
    outf = os.path.join(dirname, prefix + '_stringtie.gtf')
    command = 'stringtie -G %s -p %s -o %s --rf %s' % (gtf, thread, outf, bam)
    run_command(command, 'Error in StringTie!')


if __name__ == '__main__':
    from docopt import docopt
    assemble(docopt(__doc__, version=__version__))
