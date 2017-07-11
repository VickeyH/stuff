#!/usr/bin/env python3

'''
Usage: dbloci.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -l=SIZE                        Size of loci flanking region.
                                   [default: 200]
    --filter=FILTER                Filter setting (rpm/height).
    --rpm=RPM                      Filter of RPM. [default: 1]
    --height=HEIGHT                Filter of cluster height. [default: 2]
'''

import os.path
import numpy as np
from collections import defaultdict
from seqlib.path import check_dir
from seqlib.ngs import check_bed

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def dbloci(options):
    '''
    Fetch bidirectionally transcribed loci
    '''
    # parse options
    folder = check_dir(options['<rampagedir>'])
    size = int(options['-l'])
    filter_flag = options['--filter']
    if filter_flag == 'rpm':
        filter = float(options['--rpm'])
    elif filter_flag == 'height':
        filter = int(options['--height'])
    # fetch rampage pairs
    peak_f = os.path.join(folder, 'rampage_peaks.txt')
    peak_bed = check_bed(peak_f)
    up_cluster = {}
    down_cluster = {}
    pairs = defaultdict(list)
    with open(peak_f, 'r') as peak:
        for ucluster in peak:  # parse upstream cluster
            chrom, _, _, _, _, ustrand, upos = ucluster.split()[:7]
            if ustrand == '+':  # upstream cluster should be minus
                continue
            if filter_cluster(ucluster, filter_flag, filter):
                continue
            uheight = int(ucluster.split()[7])
            up_id = '\t'.join([chrom, upos, ustrand])
            up_cluster[up_id] = uheight
            # parse downstream cluster
            start = int(upos)
            end = start + size * 2
            for dcluster in peak_bed.fetch(chrom, start, end):
                dstrand, dpos = dcluster.split()[5:7]
                if dstrand == '-':  # downstream cluster should be plus
                    continue
                if filter_cluster(dcluster, filter_flag, filter):
                    continue
                dheight = int(dcluster.split()[7])
                down_id = '\t'.join([chrom, dpos, dstrand])
                down_cluster[down_id] = dheight
                # construct pairs
                pairs[up_id].append(down_id)
                pairs[down_id].append(up_id)
    # output enhancers
    outf = os.path.join(folder, 'enhancers.txt')
    with open(outf, 'w') as out:
        for pair_set in fetch_pair(pairs):
            up_site, down_site = 0, 0
            up_height, down_height = 0, 0
            for site_id in pair_set:
                chrom, site, strand = site_id.split()
                if strand == '-':  # upstream
                    height = up_cluster[site_id]
                    if height > up_height:
                        up_site = int(site)
                        up_height = height
                else:  # downstream
                    height = down_cluster[site_id]
                    if height > down_height:
                        down_site = int(site)
                        down_height = height
            middle_site = int((up_site + down_site) / 2)
            forward_plus = cal_density(folder, chrom, middle_site,
                                       middle_site + size, 'plus')
            forward_minus = cal_density(folder, chrom, middle_site,
                                        middle_site + size, 'minus')
            reverse_plus = cal_density(folder, chrom, middle_site - size,
                                       middle_site, 'plus')
            reverse_minus = cal_density(folder, chrom, middle_site - size,
                                        middle_site, 'minus')
            if forward_minus >= forward_plus or reverse_plus >= reverse_minus:
                continue
            else:
                forward = forward_plus
                reverse = reverse_minus
            forward_dis = fetch_dis(folder, chrom, middle_site,
                                    middle_site + size, '+')
            reverse_dis = fetch_dis(folder, chrom, middle_site - size,
                                    middle_site, '-')
            fold = (forward - reverse) * 1.0 / (forward + reverse)
            start = middle_site - size
            end = middle_site + size
            out_format = '%s\t%d\t%d\tenhancer\t0\t+\t%d\t%d\t%d\t%d\t%d\t%f\n'
            out.write(out_format % (chrom, start, end, middle_site,
                                    reverse_dis, forward_dis, reverse, forward,
                                    fold))


def fetch_dis(folder, chrom, start, end, strand):
    link_bed = check_bed(os.path.join(folder, 'rampage_link.bed'))
    dis_lst = []
    for link in link_bed.fetch(chrom, start, end):
        _, link_start, link_end, _, _, link_strand = link.split()[:6]
        if link_strand != strand:
            continue
        link_start = int(link_start)
        link_end = int(link_end)
        if strand == '+':
            if link_start < start or link_start > end:
                continue
        else:
            if link_end < start or link_end > end:
                continue
        dis_lst.append(link_end - link_start)
    return np.percentile(dis_lst, 75)


def cal_density(folder, chrom, start, end, flag):
    if flag == 'plus':  # plus
        tag_bed = check_bed(os.path.join(folder, 'rampage_plus_5end.bed'))
    else:  # minus
        tag_bed = check_bed(os.path.join(folder, 'rampage_minus_5end.bed'))
    total = 0
    for tag in tag_bed.fetch(chrom, start, end):
        total += 1
    return total


def fetch_pair(pairs):
    removed_item = set()
    for init_item in pairs:
        if init_item in removed_item:
            continue
        removed_item.add(init_item)
        pair_set = set([init_item] + pairs[init_item])
        checked_item = pairs[init_item][:]
        while checked_item:
            item = checked_item.pop()
            if item in removed_item:
                continue
            removed_item.add(item)
            pair_set.update(pairs[item])
            checked_item.extend(set(pairs[item]).difference(removed_item))
        else:
            yield pair_set


def filter_cluster(cluster, filter_flag, filter):
    if filter_flag == 'rpm':  # if filter is rpm
        rpm = float(cluster.rstrip().split()[-1])
        if rpm < filter:  # not enough rpm
            return True
    elif filter_flag == 'height':  # if filter is height
        height = int(cluster.split()[7])
        if height < filter:  # not enough height
            return True
    return False


if __name__ == '__main__':
    from docopt import docopt
    dbloci(docopt(__doc__, version=__version__))
