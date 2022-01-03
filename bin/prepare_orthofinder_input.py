#!/usr/bin/env python3

import sys
import re
import os.path
import gzip
from Bio import SeqIO, SeqRecord

species2id = {}
sequence2id = {}
blast_regex = re.compile("^(.*)_vs_(.*).txt.gz$")

with open("SpeciesIDs.txt", 'w') as species:
  with open(sys.argv[1]) as file:
    for i, line in enumerate(file):
      path = line.strip()
      name = os.path.basename(path)
      species2id[name] = i
      species.write("{0}: {1}\n".format(i, path))
      sequence2id[i] = {}
      with open("Species" + str(i) + ".fa", 'w') as fasta:
        for j, record in enumerate(SeqIO.parse(path, "fasta")):
          sequence2id[i][record.id] = str(i) + "_" + str(j)
          record.id = str(i) + "_" + str(j)
          record.description = ''
          fasta.write(record.format('fasta'))

with open("SequenceIDs.txt", 'w') as output:
  for path, u in sequence2id.items():
    for id, label in u.items():
      output.write("{0}: {1}\n".format(label, id))

with open(sys.argv[2]) as file:
  for line in file:
    path = line.strip()
    match = blast_regex.match(path)
    if match == None:
      raise ValueError('Invalid filename ' + path)
    left_id = species2id[match.group(1)]
    right_id = species2id[match.group(2)]
    with gzip.open('Blast{}_{}.txt.gz'.format(left_id, right_id), 'w') as output:
      with gzip.open(path) as input:
        for line in input:
          u = line.strip().decode().split('\t')
          u[0] = sequence2id[left_id][u[0]]
          u[1] = sequence2id[right_id][u[1]]
          output.write(('\t'.join(u) + '\n').encode())
