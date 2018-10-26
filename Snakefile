#! /usr/bin/env python
# -*- coding: utf-8 -*-
from future import standard_library
standard_library.install_aliases()

from os import path

import snakemake
from snakemake.exceptions import RuleException
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

""" Functions """
def dict_args(d):
    """ Convert a dictionary to command line arguments"""
    ret = ''
    for k,v in d.items():
        if len(k) == 1:
            ret += '-{} {} '.format(k,v)
        else:
            ret += '--{} {} '.format(k,v)
    return ret

def read_value(fn):
    """ Read a value from text file.
        
        The file should contain a single word indicating the value. For example, the 
        library layout of a run could be stored in a file named 'LAYOUT.txt', and the 
        value should be either 'PAIRED' or 'SINGLE'.
    """
    with open(fn, 'rU') as fh:
        return fh.read().strip()

configfile: "config.yaml"

wildcard_constraints:
    sampid="GSM\d+",
    sraid="[SED]RR\d+",
    sraproj="[SED]RP\d+",



''' Constants '''
SAMP_TABLE = [l.strip('\n').split('\t') for l in open(config['sample_table'], 'rU')][1:]
SAMPLES = [r[0] for r in SAMP_TABLE] # [l.strip() for l in open('cur_samples.txt', 'rU')]
RUNS = {r[0]:r[1].split(',') for r in SAMP_TABLE}

# print(RUNS)

localrules: all, complete_sample

rule all:
    input:
        expand("samples/{s}/completed.txt", s=SAMPLES)

rule complete_sample:
    input:
        "samples/{sampid}/bt2_multi.unsorted.bam",
        "samples/{sampid}/bt2_multi.summary.txt",       
        "samples/{sampid}/inform-telescope_report.tsv",
        "samples/{sampid}/RepEnrich_fraction_counts.txt",
        "samples/{sampid}/bt_repenrich.summary.txt",
        "samples/{sampid}/tetx_count.cntTable"
    output:
        touch("samples/{sampid}/completed.txt")

''' References '''
include: "Snakefile.references"

''' Run-level download and preprocessing '''
include: "Snakefile.rundata"

''' Sample-level processing '''
include: "Snakefile.sampledata"

''' Alignment (run-level)'''
include: "Snakefile.alignment"

''' Telescope '''
include: "Snakefile.telescope"

''' RepEnrich '''
include: 'Snakefile.repenrich'

''' RepEnrich '''
include: 'Snakefile.tetranscripts'
