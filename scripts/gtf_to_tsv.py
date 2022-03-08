#! /usr/bin/env python
import sys
import argparse

import re

def gtf_to_tsv(infile, outfile):
    colnames = ['locus', 'chrom', 'start', 'end', 'strand', 'score',
                'category', 'gene_id', 'intModel', 'locid', 'model_cov', 'model_pct',
                'transcript_id']
    lines = (l.strip('\n').split('\t') for l in open(infile, 'rU'))
    alld = []
    for l in lines:
        if l[2] == 'gene':
            alld.append(
                {'chrom': l[0], 'start': l[3], 'end': l[4], 
                 'strand': l[6], 'score': l[5]}
            )
            alld[-1].update({k:v for k,v in re.findall('(\S+)\s+"([\s\S]*?)";', l[8])})
    
    with open(outfile, 'w') as outh:
        print('\t'.join(colnames), file=outh)
        for d in alld:
            print('\t'.join(d[k] if k in d else '.' for k in colnames), file=outh)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Output telebuilder GTF as a TSV")
    parser.add_argument('infile', 
                        help='Input GTF file.')
    parser.add_argument('outfile', 
                        help='Output TSV file.')
    args = parser.parse_args()
    gtf_to_tsv(args.infile, args.outfile)
