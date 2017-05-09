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
from collections import defaultdict
from seqlib.path import which
from seqlib.helper import run_command

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def assemble(options):
    '''
    Assemble RNA-seq with StringTie
    '''
    # parse options
    if not which('stringtie'):
        sys.exit('Error: No StringTie installed!')
    bam = options['<bam>']
    gtf = options['--gtf']
    thread = options('--thread')
    # run StringTie
    out_gtf = run_stringtie(bam, gtf, thread)
    convert_gtf(out_gtf)


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
    return outf


def convert_gtf(out_gtf):
    # set variables
    prefix = os.path.splitext(out_gtf)[0]
    annotated_gene = defaultdict(list)
    novel_isoform = {}
    annotated_exon = defaultdict(dict)
    novel_exon = defaultdict(set)
    # parse GTF
    with open(out_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom, _, etype, start, end, _, strand, _, info = line.split('\t')
            info = {m: n.strip('"') for m, n in map(lambda x: x.split(),
                                                    info.split(';')[:-1])}
            if etype == 'transcript':  # transcript
                out_info = '\t'.join(['%s', chrom, strand, start, end, start,
                                      start, '1', start + ',', end + ',',
                                      info['FPKM'], info['TPM']]) + '\n'
                if 'ref_gene_name' in info:
                    gene_name = info['ref_gene_name']
                    annotated_gene[gene_name].append(out_info % gene_name)
                else:
                    iso_id = info['transcript_id']
                    novel_isoform[iso_id] = out_info
            else:  # exon
                if 'ref_gene_name' in info:
                    gene_id = info['gene_id']
                    gene_name = info['ref_gene_name']
                    if gene_name not in annotated_exon[gene_id]:
                        annotated_exon[gene_id][gene_name] = set([start, end])
                    else:
                        annotated_exon[gene_id][gene_name].update([start, end])
                else:
                    iso_id = info['transcript_id']
                    novel_exon[iso_id].update([start, end])
    # assign gene name
    for iso in novel_isoform:
        iso_info = novel_isoform[iso]
        iso_exons = novel_exon[iso]
        iso_exon_num = len(iso_exons)
        gene_id = iso.rsplit('.', 1)[0]
        gene_exon_set = annotated_exon[gene_id]
        for gene_name in gene_exon_set:
            gene_exons = gene_exon_set[gene_name]
            overlapped_num = len(iso_exons.intersection(gene_exons))
            if overlapped_num >= iso_exon_num * 0.5:  # have gene_name
                annotated_gene[gene_name].append(iso_info % gene_name)
                break
        else:  # no gene_name
            annotated_gene[iso].append(iso_info % iso)
    # output GenePred
    outf = prefix + '.txt'
    with open(outf, 'w') as out:
        for gene_name in annotated_gene:
            for iso_info in annotated_gene[gene_name]:
                out.write(iso_info)


if __name__ == '__main__':
    from docopt import docopt
    assemble(docopt(__doc__, version=__version__))
