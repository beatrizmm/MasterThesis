#!/usr/bin/env python

import sys
import re

ref_file = sys.argv[1]    #Reference transcript lengths
read_file = sys.argv[2]   #Read lengths
pair_file = sys.argv[3]    #Best match pairs list (already filtered by 25%) obtained with Evaluation_nucleotide_level.py

reflen = {}     #Dictionary with ID as keys and total length as value
readlen = {}


def parse_file(lengths_file): 
  idlen = {}
  with open(lengths_file, "r") as len_list: 
    for i in len_list: 
      line = i.split()
      idlen[line[1]] = line[0]
  len_list.close()
  return (idlen)

reflen = parse_file(ref_file)
readlen = parse_file(read_file)

with open (pair_file, "r") as pairs:
  for i in pairs: 
    line = i.split()
    proportion = int(readlen[line[0]])/(int(reflen[line[1]])*1.0)
    printout = [line[0], line[1], format(proportion, '.4f')]
    print ('\t'.join(printout))
pairs.close()

