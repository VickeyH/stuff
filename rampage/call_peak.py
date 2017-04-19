#!/usr/bin/env python3

'''
Usage: call_peak.py [options] -g  -5 PEAK5 -3 PEAK3 <rampagebed>
Options:
    -h --help                      Show help message.
    --version                      Show version.
    -5 PEAK5                       F-seq peak file of rampage 5end.
    -3 PEAK3                       F-seq peak file of rampage 3read.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output prefix. [default: rampage]
'''

import math
import numpy as np
import scipy.stats

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def call_peak(options):
    '''
    Call rampage peaks
    '''
    pass
    # parse options
    # parse gene annotation
    # align candidate peak
    # filter peaks


def wald_test(value, samples):
    '''
    Test whether value belongs to samples
    '''
    se = math.sqrt(np.var(samples) / samples.size)
    w = (value - np.mean(samples)) / se
    p = scipy.stats.norm.sf(w)
    return w, p


if __name__ == '__main__':
    from docopt import docopt
    call_peak(docopt(__doc__, version=__version__))
