#! /usr/bin/env python
import sys
import argparse

import pysam

def fetch_bundle(samfile, **kwargs):
    samiter = samfile.fetch(**kwargs)
    bundle = [ next(samiter) ]
    for aln in samiter:
        if aln.query_name == bundle[0].query_name:
            bundle.append(aln)
        else:
            yield bundle
            bundle = [aln]
    yield bundle


def make_unique_bam(inname, outname, layout):
    inbam = pysam.AlignmentFile(inname, 'rb')

    with pysam.AlignmentFile(outname, 'wb', template=inbam) as outbam:
        if layout == 'PAIRED':
            for bundle in fetch_bundle(inbam, until_eof=True):
                if len(bundle) == 2:
                    if bundle[0].is_unmapped and bundle[1].is_unmapped:
                        continue
                    outbam.write(bundle[0])
                    outbam.write(bundle[1])
        else:
            for bundle in fetch_bundle(inbam, until_eof=True):
                if len(bundle) == 1:
                    if bundle[0].is_unmapped:
                        continue
                    outbam.write(bundle[0])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Output uniquely mapping read from a multimapped alignment")
    parser.add_argument('--layout', choices=['PAIRED','SINGLE'], default='PAIRED',
                        help='Library layout')
    parser.add_argument('infile', 
                        help='Input BAM file (Multimapped)')
    parser.add_argument('outfile', 
                        help='Output BAM file (Uniquely mapped)')
    args = parser.parse_args()
    make_unique_bam(args.infile, args.outfile, args.layout)
