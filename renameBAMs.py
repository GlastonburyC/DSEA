#!/usr/bin/env python

import os,subprocess
import glob

mapping = {}

# load in the mapping file [SubjectID    BAM_id]
with open('/gpfs/home/DTR/Expression/EUROBATS/DataRelease/keySamples.F') as mapping_file:
    for line in (line.strip().split() for line in mapping_file):
        mapping[line[2]] = line[0]

suffix = '_sorted.bam.map10.f3.NM2.sorted'

filename = glob.glob('*.bam')
for i in filename:
    root, extension = os.path.splitext(i)
    stripped_root = root[:-len(suffix)]
    if stripped_root in mapping:
        os.rename(i, ''.join(mapping[stripped_root] + '.bam'))


with open('/gpfs/home/DTR/Expression/EUROBATS/Counts/qc_F_freezev1.txt') as mapping_file:
    for line in (line.strip().split() for line in mapping_file):
        if os.path.isfile(line[0]+'.bam'):
            subprocess.call(["mv",line[0]+'.bam', "bams_with_genotypes/"])
