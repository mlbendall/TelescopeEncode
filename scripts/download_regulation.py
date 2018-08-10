#! /usr/bin/env python
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
import sys

lines = [l.strip('\n').split('\t') for l in open('regulation/ENCODE_regulation.tsv', 'rU')]
header = lines[0]
lines = lines[1:]

celltypes = sorted(set(l[0] for l in lines))
for ct in celltypes:
    print('mkdir -p regulation/%s' % ct)

url_base = 'https://www.encodeproject.org/files/%s/@@download/%s.bed.gz'

for l in lines:
    ct = l[0]
    assay = l[1]
    id = l[3]
    if id == '.':
        print('no file: %s %s' % (ct, assay), file=sys.stderr)
        continue
    dest = 'regulation/%s/%s.%s.bed.gz' % (ct, assay, id)
    dest_short = 'regulation/%s/%s.bed.gz' % (ct, assay)
    url = url_base % (id, id)
    print('[[ ! -e "%s" ]] && wget -O %s %s' % (dest, dest, url))
    print('rm -f %s && ln -s %s.%s.bed.gz %s' % (dest_short, assay, id, dest_short))
