#!/usr/bin/env python

import sys
import re
import os

ref_annot = sys.argv[1]    #Reference annotation (filtered only for 1 type of feature: Exon, CDS...)
pred_annot = sys.argv[2]   #Predicted annotation (filtered only for 1 type of feature: Exon, CDS...)
ref_table = sys.argv[3]    #List of transcripts and the gene they correspond to. This is specially for GENCODE since they do not have the IDs with GENE.T1. Structure is transID tab GeneID
trans_fasta = sys.argv[4]    #File obtained in the step05 of the annotation that contains the fasta sequences for the transcripts (LABEL.transcripts.fa)
label = sys.argv[5]         

#---------------FUNCTIONS------------------------------------------------------------------------------------------------------------------------------------

def parse_gff(gff_file):
  dict_tID = {}   
  with open (gff_file, "r") as gff:
    for i in gff: 
      line = i.split()
      ptID = re.search(r'(?<=Parent=)[^;]+', line[8]).group()
      length = int(line[4]) - int(line[3]) +1
      if ptID in dict_tID: 
        dict_tID[ptID] = dict_tID[ptID] + length  
      else: 
        dict_tID[ptID] = length
  gff.close()
  return (dict_tID)
  #This creates a dictionary with tID as key and a summatory of the exon bases of each transcript.

def parse_lookuptbl(tbl_file): 
  dict_trans_gene = {}
  with open (tbl_file, "r") as tbl: 
    for i in tbl: 
      pair = i.split()
      dict_trans_gene[pair[0]] = pair[1]
  tbl.close()
  return (dict_trans_gene)
  #This creates a dictionary with transcriptID as key and its geneID as the value

def parse_fastalen(fastalen_tbl): 
  dict_len_tID = {}
  with open (fastalen_tbl, "r") as tbl: 
    for i in tbl: 
      entry = i.split()
      dict_len_tID[entry[1]] = int(entry[0])
  tbl.close()
  return (dict_len_tID)

def matching(dict_refID, dict_ref, dict_pred, current_tID, list_genetracker): 
  best_overlap = 0
  best_match = ""
  multiple_best = []
  for i in list_genetracker: 
    if dict_refID[i] > best_overlap:
      best_overlap =  dict_refID[i]
      best_match = i
      multiple_best = [i]
# This way I get the highest overlap, but I will go twice through the list in order to check if there are multiple hits with the same score

  for x in list_genetracker: 
    if (best_overlap == dict_refID[x] and x != best_match): 
      multiple_best.append(x)

  if len(multiple_best) > 1: 
    len_proportion = []
    for m in multiple_best:
      proportion = dict_ref[m]/(dict_pred[current_tID] * 1.0) #Forcing it to be float
      len_proportion.append([m, proportion])
      list_proportions = []
    for p in len_proportion: 
      list_proportions.append(p[1])
    best_proportion = min(list_proportions, key=lambda x:abs(x-1))     #Chooses the value closest to 1
    for c in len_proportion:
      if best_proportion == c[1]: 
        predID_match = [current_tID, c[0], best_overlap]
#If multiple hits have the same score, we make the decision based on their length being as close to the ref as possible.

  elif len(multiple_best) == 1: 
    predID_match = [current_tID, best_match, best_overlap]

  return (predID_match)

#------------------------------------------------------------------------------------------------------------------------------------------------------

dict_ref = parse_gff(ref_annot)
dict_pred = parse_gff(pred_annot)
dict_reftbl = parse_lookuptbl(ref_table)
BT_output = label+".BT.intersect.out"
cmd = "intersectBed -wo -s -a "+pred_annot+" -b "+ref_annot+" > "+BT_output
os.system(cmd)
dict_predID = {}
dict_refID = {}

with open (BT_output, "r") as BTout: 
  for i in BTout: 
    line = i.split()
    predID = re.search(r'(?<=Parent=)[^;]+', line[8]).group()
    refID = re.search(r'(?<=Parent=)[^;]+', line[17]).group()
    overlap = int(line[18])
    if predID not in dict_predID:
      dict_predID[predID] = []
      #print predID  ######

      #Since it is possible for a reference transcript to be a hit of multiple prediction transcript, I empty and restart dict_refID for each predID. 
      #This is done assuming that the bedtools output is ordered by the gff file given as -a.

      if dict_refID:      #Avoids making the match for the best isoform in the first line of the file
        counter = 1
        for g in gene_tracker: 
          predID_match = matching(dict_refID, dict_ref, dict_pred, current_tID, gene_tracker[g])
          if len(gene_tracker) == 1: 
            dict_predID[predID_match[0]] = [predID_match[1], predID_match[2]]
          elif len(gene_tracker) > 1: 
            new_predID = predID_match[0]+'*'+str(counter)
            dict_predID[new_predID] = [predID_match[1], predID_match[2]]
            counter = counter +1
            
#In the case that a same pred transcript has multiple ref hits (because they belong to different genes), they cannot be added as they are to the dict because their ID will not be unique anymore.
#What I proposed is changing their ID by adding a * and a number. 

      dict_refID = {}  
      current_tID = predID
      dict_refID[refID] = overlap
      gene_tracker = {}
      gene_tracker[dict_reftbl[refID]]=[refID]

    else:     
      if refID in dict_refID:
        dict_refID[refID] = dict_refID[refID] + overlap
        #print current_tID, refID, dict_refID[refID] #####
      else: 
        dict_refID[refID] = overlap
        #print current_tID, refID, dict_refID[refID] #####
        if dict_reftbl[refID] in gene_tracker:
          gene_tracker[dict_reftbl[refID]].append(refID)
          #print dict_reftbl[refID], gene_tracker[dict_reftbl[refID]]
        else: 
          gene_tracker[dict_reftbl[refID]] = [refID]
          #print dict_reftbl[refID], gene_tracker[dict_reftbl[refID]]
            
BTout.close()
    
#The last transcript of the file needs to be handled outside the loop (because it won't access the best matching step). Bad practices, I know, I should define a function.

counter = 1
for g in gene_tracker: 
  predID_match = matching(dict_refID, dict_ref, dict_pred, current_tID, gene_tracker[g])
  if len(gene_tracker) == 1: 
    dict_predID[predID_match[0]] = [predID_match[1], predID_match[2]]
  elif len(gene_tracker) > 1: 
    new_predID = predID_match[0]+'*'+counter
    dict_predID[new_predID] = [predID_match[1], predID_match[2]]
    counter = counter +1

#Now we need to calculate SP, SN and accuracy for each pair of transcripts

eval_output = label+".evaluation_nt_level.tbl"
final_output = label+".evaluation_nt_level.filtered.tbl"   #This is the filtered output that we will obtain at the end but I amb already writing the ones without match

with open (eval_output, "w") as EVALout:
  with open (final_output, "w") as EVALfilteredout:
    #pr_title = ["Predicted transcript", "Reference transcript", "SP", "SN", "AC\n"]
    #EVALout.write('\t'.join(pr_title))
    for t in dict_pred: 
      if t in dict_predID: 
        if not dict_predID[t]:    #If the list for best match is empty means that there were multiple hits in different genes, so all the info will be with tID*counter. 
          for b in dict_predID: 
            if (t+"*") in b: 
              sp = dict_predID[b][1] / (dict_pred[t] * 1.0)
              sn = dict_predID[b][1] / (dict_ref[dict_predID[b][0]] * 1.0)
              ac = (sp + sn)/2
              pr_out = [t, dict_predID[b][0], format(sp, '.4f'), format(sn, '.4f'), format(ac, '.4f')+"\n"]    #Print tID with possible duplicated entries w/ multiple gencode hits
              #pr_out = [b, dict_predID[b][0], format(sp, '.4f'), format(sn, '.4f'), format(ac, '.4f')+"\n"]   #Print the tID to make it unique
              EVALout.write('\t'.join(pr_out))
            
        else: 
          sp = dict_predID[t][1] / (dict_pred[t] * 1.0)
          sn = dict_predID[t][1] / (dict_ref[dict_predID[t][0]] * 1.0)
          ac = (sp + sn)/2
          pr_out = [t, dict_predID[t][0], format(sp, '.4f'), format(sn, '.4f'), format(ac, '.4f')+"\n"]
          EVALout.write('\t'.join(pr_out))

      else: 
        pr_out = [t, '-', "0", "0", "0\n"]
        EVALout.write('\t'.join(pr_out))
        EVALfilteredout.write('\t'.join(pr_out))

EVALout.close()
EVALfilteredout.close()

#Then we filter the results to get only the best hit from the reference, based on best accuracy and from the ones with multiple best, the longest one. 

#First we filter the output
output_filtered = label+".evaluation_nt_level.sorted.tbl"
cmd = "grep -v '-' "+eval_output+" | sort -k2,2 -k5,5nr > "+output_filtered
os.system(cmd)

#Second we obtain the lengths of the transcripts (we do this with external file because if we are comparing CDS we want to take into account the UTRs to decide the best match between same AC)
output_lengths = label+".fastalength.tbl"
cmd = "fastalength "+trans_fasta+" > "+output_lengths
os.system(cmd)
dict_predID_len = parse_fastalen(output_lengths)


dict_best_ref_hit = {}
with open (output_filtered, "r") as TBLout: 
  for k in TBLout: 
    line = k.split()
    refID = line[1]
    if refID in dict_best_ref_hit: 
      if float(line[4]) == dict_best_ref_hit[refID][3]: 
        if dict_predID_len[line[0]] > dict_predID_len[dict_best_ref_hit[refID][0]]: 
          dict_best_ref_hit[refID] = [line[0], float(line[2]), float(line[3]), float(line[4])]
    else: 
      dict_best_ref_hit[refID] = [line[0], float(line[2]), float(line[3]), float(line[4])]
TBLout.close()

with open (final_output, "a") as EVALfilteredout:
  for p in dict_best_ref_hit: 
    pr_line = [dict_best_ref_hit[p][0], p, str(dict_best_ref_hit[p][1]), str(dict_best_ref_hit[p][2]), str(dict_best_ref_hit[p][3])+"\n"]
    EVALfilteredout.write('\t'.join(pr_line))

EVALfilteredout.close()








