#!/usr/bin/env python

import sys
import re
import os

ref_annot = sys.argv[1]    #Reference annotation (filtered only for one type of feature: CDS, exon...)
pred_annot = sys.argv[2]   #Predicted annotation (filtered only for one type of feature: CDS, exon...)
eval_from_nt_output = sys.argv[3]   #This is the output file obtained from running the Evaluation_nucleotide_level, the .tbl file
label = sys.argv[4]         

#---------------FUNCTIONS-------------------------

def parse_gff(gff_file):
  dict_tID_coord = {}  
  with open (gff_file, "r") as gff:
    for i in gff: 
      line = i.split()
      ptID = re.search(r'(?<=Parent=)[^;]+', line[8]).group()
      if ptID in dict_tID_coord: 
        dict_tID_coord[ptID].append([line[3], line[4]])
      else: 
        dict_tID_coord[ptID]=[]
        dict_tID_coord[ptID].append([line[3], line[4]])
  gff.close()
  return (dict_tID_coord)
  #This creates a a dictionary with the transcript ID as key and a list of lists as values, where each list is an exon with its start and end positions. 

def parse_tbl(tbl_file):      
  dict_matches = {}
  dict_nomatch = {}
  with open (tbl_file, "r") as matches: 
    for i in matches: 
      line = i.split()
      if line[1] == "-": 
        dict_nomatch[line[0]] = line[1]
      else: 
        dict_matches[line[1]] = line[0]
  matches.close()
  return (dict_matches, dict_nomatch)
#Taking into account that this tbl file has unique hits for the reference annotation (coming from the eval_nt_level_filtered output), since the predicted tID are not unique. 
#I create two dictionaries, one with the matches where refID are the keys and predID are the values, and the other one with predID as keys since there are no hits for them. 


#-------------------------------------------------

dict_ref_coord = parse_gff(ref_annot)
dict_pred_coord = parse_gff(pred_annot)
(dict_matches, dict_nomatches) = parse_tbl(eval_from_nt_output)

#dict_starts = {}
#dict_ends = {}
dict_perfect = {}

for i in dict_matches:
  dict_perfect[i] = 0
  for e in dict_pred_coord[dict_matches[i]]: 
    ####print e
    for r in dict_ref_coord[i]:
      if (e[0] == r[0] and e[1] == r[1]):           #Only considering perfect matches, when start coordinate and end coordinate match. 
        dict_perfect[i] = dict_perfect[i] +1


#Now we need to calculate SP, SN and accuracy

eval_output = label+".evaluation_exon_level.tbl"
with open (eval_output, "w") as EVALout:
  #pr_title = ["Predicted transcript", "Reference transcript", "SP", "SN", "AC\n"]
  #EVALout.write('\t'.join(pr_title))
  for t in dict_matches: 
    pred_exons = len(dict_pred_coord[dict_matches[t]])
    ref_exons = len(dict_ref_coord[t])
    sp = dict_perfect[t] / (pred_exons * 1.0)
    sn = dict_perfect[t] / (ref_exons * 1.0)
    ac = (sp + sn)/2
    pr_out = [dict_matches[t], t, format(sp, '.4f'), format(sn, '.4f'), format(ac, '.4f')+"\n"]
    EVALout.write('\t'.join(pr_out))

  for m in dict_nomatches: 
    pr_out = [m, dict_nomatches[m], "0", "0", "0\n"]
    EVALout.write('\t'.join(pr_out))

EVALout.close()



#for i in dict_ref: 
#  print i








