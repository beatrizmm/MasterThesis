#!/usr/bin/env python

import sys
import re

ref_annot = sys.argv[1]    #Reference annotation
pred_annot = sys.argv[2]   #Predicted annotation
eval_from_nuc_output = sys.argv[3]   #This is the output file obtained from running the Evaluation_nucleotide_level.py, the filtered.tbl file. IMPORTANT! Better using the comparison between CDS. 
label = sys.argv[4]         

#--------------------------------FUNCTIONS-----------------------------------------------------------------------------------------------------------------------------

def parse_gff(gff_file):
  dict_pos_trans = {}  
  dict_neg_trans = {}
  dict_exon = {}         #This will be a list of lists with all the exons in the transcript, start and end coordinate
  dict_CDS = {}          #This will be a unique list, of the smallest coordinate of CDS and largest coordinate of CDS
  with open (gff_file, "r") as gff:
    for i in gff: 
      line = i.split()
      if line[2] == "transcript":
        tID =  re.search(r'(?<=ID=)[^;]+', line[8]).group()
        if line[6] == "+": 
          dict_pos_trans[tID] = []
        elif line[6] == "-": 
          dict_neg_trans[tID] = []

      elif line[2] == "exon": 
        if tID in dict_exon: 
          dict_exon[tID].append([int(line[3]), int(line[4])])
        else: 
          dict_exon[tID] = []
          dict_exon[tID].append([int(line[3]), int(line[4])])

      elif line[2] == "CDS": 
        if tID in dict_CDS: 
          if int(line[3]) < dict_CDS[tID][0]:
            dict_CDS[tID][0] = int(line[3])
          if int(line[4]) > dict_CDS[tID][1]:
            dict_CDS[tID][1] = int(line[4])
        else: 
          dict_CDS[tID] = [int(line[3]), int(line[4])]
  gff.close()

  dict_3UTR = {}
  dict_5UTR = {}

  #Transcripts in the positive strand
  for d in dict_pos_trans: 
    dict_5UTR[d] = []
    dict_3UTR[d] = []
    for t in dict_exon[d]: 
      #5' UTR evaluation
      if t[0] < dict_CDS[d][0]: 
        if t[1] < dict_CDS[d][0]: 
          dict_5UTR[d].append([t[0],t[1]])
        else: 
          dict_5UTR[d].append([t[0], dict_CDS[d][0]])

      #3' UTR evaluation
      if t[1] > dict_CDS[d][1]: 
        if t[0] > dict_CDS[d][1]: 
          dict_3UTR[d].append([t[0], t[1]])
        else: 
          dict_3UTR[d].append([dict_CDS[d][1], t[1]])

  #Transcripts in the negative strand
  for n in dict_neg_trans: 
    dict_5UTR[n] = []
    dict_3UTR[n] = []
    for a in dict_exon[n]: 
      #5' UTR evaluation
      if a[1] > dict_CDS[n][1]:
        if a[0] > dict_CDS[n][1]: 
          dict_5UTR[n].append([a[0],a[1]])
        else: 
          dict_5UTR[n].append([dict_CDS[n][1], a[1]])

      #3' UTR evaluation
      if a[0] < dict_CDS[n][0]: 
        if a[1] < dict_CDS[n][0]: 
          dict_3UTR[n].append([a[0], a[1]])
        else: 
          dict_3UTR[n].append([a[0], dict_CDS[n][0]])

  return (dict_5UTR, dict_3UTR)


def parse_tbl(tbl_file):                #Since tID are not unique, the key of the dict will be the reference ID and the values the predID
  dict_matches = {}
  with open (tbl_file, "r") as matches: 
    for i in matches: 
      line = i.split()
      if line[1] != "-":                #Only compare UTRs of the ones that actually have a match, DUH
        dict_matches[line[1]] = line[0]
  matches.close()
  return (dict_matches)


def UTR_len_dif(dict_pred, dict_ref, IDpred, IDref): 
  len_UTR_pred = 0
  len_UTR_ref = 0
  dif_UTR = 0
  if dict_pred[IDpred]:
    sum_pred = 0
    for i in dict_pred[IDpred]: 
      sum_pred = i[1] - i[0] +1
      len_UTR_pred = len_UTR_pred + sum_pred
  if dict_ref[IDref]:
    sum_ref = 0
    for j in dict_ref[IDref]: 
      sum_ref = j[1] - j[0] +1
      len_UTR_ref = len_UTR_ref + sum_ref

  dif_UTR = len_UTR_pred - len_UTR_ref
  return (dif_UTR)


#--------------------------------------------------------------------------------------------------------------------------------------------------------

(dict_5UTR_ref, dict_3UTR_ref) = parse_gff(ref_annot)
(dict_5UTR_pred, dict_3UTR_pred) = parse_gff(pred_annot)
dict_matches = parse_tbl(eval_from_nuc_output)

#Evaluate length difference between the pairs of transcripts

eval_output = label+".evaluation_UTR_diff.tbl"
with open (eval_output, "w") as EVALout:
  #pr_title = ["Predicted transcript", "Reference transcript", "5'UTR", "3'UTR", "Combined UTR\n"]
  #EVALout.write('\t'.join(pr_title))
  for m in dict_matches:
    #For 5'UTR
    dif_5UTR = UTR_len_dif(dict_5UTR_pred, dict_5UTR_ref, dict_matches[m], m)

    #For 3'UTR
    dif_3UTR = UTR_len_dif(dict_3UTR_pred, dict_3UTR_ref, dict_matches[m], m)

    #Combined
    combined = dif_5UTR + dif_3UTR

    pr_out = [dict_matches[m], m, str(dif_5UTR), str(dif_3UTR), str(combined)+"\n"]
    EVALout.write('\t'.join(pr_out))

EVALout.close()















