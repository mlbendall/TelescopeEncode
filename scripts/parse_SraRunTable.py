#! /usr/bin/env python
from future import standard_library
standard_library.install_aliases()

from collections import defaultdict, Counter

lines = (l.strip('\n').split('\t') for l in open('metadata/SraRunTable.txt', 'rU'))

# Get header and add fields for 'subcellular_fraction' and 'sequencing_center'
header = next(lines)
header = header + ['subcellular_fraction', 'sequencing_center']

newlines = []
for l in lines:
    lib_name = l[header.index('Library_Name')].split()[1]
    # Determine RNA extract
    rnaextract = l[header.index('rnaextract')]
    if rnaextract == '':
        assert lib_name.split('_')[0] == 'Caltech' # All Caltech are longPolyA
        l[header.index('rnaextract')] = 'longPolyA'
    # Determine subcellular fraction
    sub_frac = 'cell' if lib_name.split('_')[0] == 'Caltech' else lib_name.split('_')[-2]
    l.append(sub_frac)
    # Determine sequencing center
    if lib_name.split('_')[0] == 'Caltech':
        l.append('CALTECH')
    else:
        assert lib_name.split('_')[0][:4].upper() == 'CSHL'
        l.append('CSHL')
    newlines.append(l)

# Gather lines by sample name
by_sample = defaultdict(list)
for l in newlines:
    by_sample[l[header.index('Sample_Name')]].append(l)

# Long sample table with all fields
sample_table = []
for sampid, runs in by_sample.items():
    sample_table.append([])
    for i,f in enumerate(header):
        if all(r[i] == runs[0][i] for r in runs):
           sample_table[-1].append(runs[0][i])
        else:
            sample_table[-1].append(','.join(r[i] for r in runs))

with open('metadata/sampletable_long.tsv', 'w') as outh:
    print('\t'.join(_ for _ in header), file=outh)
    for sline in sample_table:
        print('\t'.join(_ for _ in sline), file=outh)

# Short sample table
short_table = []
short_header = ['Sample_Name', 'Run', 'source_name', 'rnaextract', 'subcellular_fraction', 
                'LibraryLayout', 'readtype', 'sequencing_center']
cols = [header.index(f) for f in short_header]
for sline in sample_table:
    short_table.append([sline[c] for c in cols])

select_types = ['H1-hESC', 'GM12878', 'K562', 
                'HeLa-S3', 'HepG2', 'HUVEC',
                'SK-N-SH', 'IMR90', 'A549', 'MCF-7', 'CD20+', 'Monocytes-CD14+', 'NHEK']

select_extract = ['longPolyA', 'total', 'longNonPolyA']

set1_table = [_ for _ in short_table if _[2] in select_types and _[4] == 'cell']
set2_table = [_ for _ in short_table if _[2] in select_types and _[4] != 'cell']
set3_table = [_ for _ in short_table if _[2] not in select_types]

with open('metadata/set1.tsv', 'w') as outh:
    print('\t'.join(short_header), file=outh)
    for l in set1_table:
        print('\t'.join(l), file=outh)

with open('metadata/set2.tsv', 'w') as outh:
    print('\t'.join(short_header), file=outh)
    for l in set2_table:
        print('\t'.join(l), file=outh)

with open('metadata/set3.tsv', 'w') as outh:
    print('\t'.join(short_header), file=outh)
    for l in set3_table:
        print('\t'.join(l), file=outh)
