#!/usr/bin/env python3

'''
Usage: assemble.py [options] -g GTF <bam>...

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -g GTF --gtf=GTF               Gene annotation GTF file.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    --dir=DIR                      Output directory. [default: ./]
    --prefix=PREFIX                Prefix for merged GTF.
                                   [default: merged_stringtie]
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
    gtf = options['--gtf']
    thread = options['--thread']
    bamf = options['<bam>']
    dir_path = options['--dir']
    if len(bamf) == 1:  # no replicates
        # run StringTie
        out_gtf = run_stringtie(bamf[0], gtf, dir_path, thread)
        # convert GTF to GenePred
        convert_gtf(out_gtf)
    else:  # have replicates
        # run StringTie
        gtf_list = []
        for f in bamf:
            gtf_list.append(run_stringtie(f, gtf, dir_path, thread))
        merged_prefix = options['--prefix']
        out_gtf = merge_stringtie(gtf_list, gtf, dir_path, merged_prefix,
                                  thread)
        # convert GTF to GenePred
        convert_gtf(out_gtf, merge=True)


def run_stringtie(bam, gtf, dir_path, thread):
    if not os.path.isfile(bam):
        sys.exit('No BAM file: %s!' % bam)
    if not os.path.isfile(gtf):
        sys.exit('No GTF file: %s' % gtf)
    fname = os.path.basename(bam)
    prefix = os.path.splitext(fname)[0]
    outf = os.path.join(dir_path, prefix + '_stringtie.gtf')
    command = 'stringtie -G %s -p %s -o %s --rf %s' % (gtf, thread, outf, bam)
    run_command(command, 'Error in StringTie!')
    return outf


def merge_stringtie(gtf_list, gtf, dir_path, prefix, thread):
    outf = os.path.join(dir_path, prefix + '.gtf')
    command = 'stringtie --merge -G %s -p %s -o %s %s' % (gtf, thread, outf,
                                                          '\t'.join(gtf_list))
    run_command(command, 'Error in StringTie --merge!')
    return outf


def convert_gtf(out_gtf, merge=False):
    # TODO: add exhausted mode
    # set variables
    prefix = os.path.splitext(out_gtf)[0]
    annotated_gene = defaultdict(list)
    novel_isoform = {}
    annotated_exon = defaultdict(dict)
    novel_exon = defaultdict(set)
    if merge:
        gene_symbol = 'gene_name'
    else:
        gene_symbol = 'ref_gene_name'
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
                                      start, '1', start + ',', end + ','])
                out_info += '\n'
                if gene_symbol in info:
                    gene_name = info[gene_symbol]
                    annotated_gene[gene_name].append(out_info % gene_name)
                else:
                    iso_id = info['transcript_id']
                    novel_isoform[iso_id] = out_info
            else:  # exon
                if gene_symbol in info:
                    gene_id = info['gene_id']
                    gene_name = info[gene_symbol]
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
        gene_id = iso.rsplit('.', 1)[0]
        gene_exon_set = annotated_exon[gene_id]
        candidate_name = ''
        max_overlap = 0
        for gene_name in gene_exon_set:
            gene_exons = gene_exon_set[gene_name]
            overlapped_num = len(iso_exons.intersection(gene_exons))
            if overlapped_num > max_overlap:
                candidate_name = gene_name
                max_overlap = overlapped_num
        if candidate_name:
            annotated_gene[candidate_name].append(iso_info % candidate_name)
        else:
            annotated_gene[gene_id].append(iso_info % gene_id)
    # output GenePred
    outf = prefix + '.txt'
    with open(outf, 'w') as out:
        for gene_name in annotated_gene:
            for iso_info in annotated_gene[gene_name]:
                out.write(iso_info)


if __name__ == '__main__':
    from docopt import docopt
    assemble(docopt(__doc__, version=__version__))
