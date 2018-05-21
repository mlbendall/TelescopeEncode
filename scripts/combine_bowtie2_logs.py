# -*- coding: utf-8 -*-
from future import standard_library
standard_library.install_aliases()

import sys
import argparse

def combine_paired(logfiles, outh):
    """
    46068964 reads; of these:
      46068964 (100.00%) were paired; of these:
        11524452 (25.02%) aligned concordantly 0 times
        25230674 (54.77%) aligned concordantly exactly 1 time
        9313838 (20.22%) aligned concordantly >1 times
        ----
        11524452 pairs aligned concordantly 0 times; of these:
          2379524 (20.65%) aligned discordantly 1 time
        ----
        9144928 pairs aligned 0 times concordantly or discordantly; of these:
          18289856 mates make up the pairs; of these:
            10385208 (56.78%) aligned 0 times
            6171063 (33.74%) aligned exactly 1 time
            1733585 (9.48%) aligned >1 times
    88.73% overall alignment rate
    """
    cols = ['run', 'reads', 'paired', 'concord_0', 'concord_1', 'concord_M',
            'x0', 'discord_1',
            'x1', 'mates', 'aln_0', 'aln_1', 'aln_M',
            'arate']
    
    newlines = []
    for f in logfiles:
        lines = [l.strip().split()[0] for l in open(f, 'rU')]
        assert lines[2] == lines[6]
        assert int(lines[2]) - int(lines[7]) == int(lines[9])
        counts = list(map(int, lines[:5] + lines[6:8] + lines[9:-1]))
        arate = lines[-1].strip('%')
        calcrate = '%.2f' % ((((counts[0] * 2) - counts[9]) / (counts[0] * 2)) * 100)
        assert calcrate == arate
        newlines.append(counts + [arate])
    
    totals = [sum(_[i] for _ in newlines) for i in range(12)]
    total_arate = '%.2f' % ((((totals[0] * 2) - totals[9]) / (totals[0] * 2)) * 100)
    
    print('\t'.join(cols), file=outh)
    print('\t'.join(['total'] + list(map(str, totals)) + [total_arate]), file=outh)
    for f,nl in zip(logfiles, newlines):
        print('\t'.join([f] + list(map(str, nl))), file=outh)


def combine_single(logfiles, outh):
    """
    21239148 reads; of these:
      21239148 (100.00%) were unpaired; of these:
        7070003 (33.29%) aligned 0 times
        9005273 (42.40%) aligned exactly 1 time
        5163872 (24.31%) aligned >1 times
    66.71% overall alignment rate
    """
    cols = ['run', 'reads', 'unpaired',
            'aln_0', 'aln_1', 'aln_M',
            'arate']
    newlines = []
    for f in logfiles:
        lines = [l.strip().split()[0] for l in open(f, 'rU')]
        counts = list(map(int, lines[:5]))
        arate = lines[-1].strip('%')
        newlines.append(counts + [arate])
    
    totals = [sum(_[i] for _ in newlines) for i in range(5)]
    total_arate = '%.2f' % ((((totals[0]) - totals[2]) / (totals[0])) * 100)
    
    print('\t'.join(cols), file=outh)
    print('\t'.join(['total'] + list(map(str, totals)) + [total_arate]), file=outh)
    for f,nl in zip(logfiles, newlines):
        print('\t'.join([f] + list(map(str, nl))), file=outh)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Output uniquely mapping read from a multimapped alignment")
    parser.add_argument('--layout', choices=['PAIRED','SINGLE'], default='PAIRED',
                        help='Library layout')
    parser.add_argument('logfiles', nargs='*',
                        help='Input Bowtie2 summary file')
    args = parser.parse_args()
    if args.layout == 'PAIRED':
        combine_paired(args.logfiles, sys.stdout)
    elif args.layout == 'SINGLE':
        combine_single(args.logfiles, sys.stdout)

