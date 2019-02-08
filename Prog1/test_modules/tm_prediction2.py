# -*- coding: utf-8 -*-
"""
tm_prediction_test.py
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""

# Your task is to *accurately* predict the primer melting points using machine 
# learning based on the sequence of the primer.

# Load the primers and their melting points.
import alignment

def ReverseCompliment(seq):
    rev_comp = ""
    for j in seq[::-1]:
        if j == "A":
            rev_comp += "T"
        if j == "T":
            rev_comp += "A"
        if j == "G":
            rev_comp += "C"
        if j == "C":
            rev_comp += "G"
    return rev_comp

def CalculatePrimerFeatures(seq):
    # modify this function to return a python list of feature values for a given sequence for Task 1

    a = seq.count("A")
    t = seq.count("T")
    c = seq.count("C")
    g = seq.count("G")

    return [a+t, c+g, len(seq)]
    
def existsmatch(primer, seq, target_score):
    (best, optloc, A) = alignment.local_align(seq, primer)
    # print("best:", best, "target score:", target_score)
    if best < target_score:
        return -1
    return optloc[0]

    
def PredictPCRProduct(sequence, fprimer, rprimer, melting_point_model):
    # modify this function for Task 2
    # all sequences should be passed in in 5' -> 3' direction

    # primers need to be 18-25 bases long
    n1 = len(fprimer)
    n2 = len(rprimer)
    if not(18 <= n1 <= 25 and 18 <= n2 <= 25):
      # print("Length test failed, n1 is %d and n2 is %d" % (n1, n2))
      return None

    # make sure both are roughly 60 degrees celsius
    forward_features = CalculatePrimerFeatures(fprimer)
    reverse_features = CalculatePrimerFeatures(rprimer)
    fT = melting_point_model.predict([forward_features])[0]
    rT = melting_point_model.predict([reverse_features])[0]
    if abs(60 - fT) > 3 or abs(60 - rT) > 3:
      # print("Melting temperature test failed, ft: %d, rt: %d" % (fT, rT))
      return None  # temperatures are more than 3 away from 60

    # make sure primers bind to sequences
    floc = existsmatch(fprimer, sequence, .9 * 10 * n1)
    if floc == -1:
      # print("Fprimer alignment test failed")
      return None

    rloc = existsmatch(ReverseCompliment(rprimer), sequence, .9*10*n2)
    if (rloc == -1):
      # print("Rprimer alignment test failed")
      return None

    endI = rloc
    startI = floc - len(fprimer)

    if abs(endI - startI) > 1000:
      # print("Binding distance test failed")
      return None

    # we're all good!
    return sequence[startI:endI]
    
def verify_PCR_exists(sequences, rprimer_list, organism_list, fprimer, melting_point_model):
  '''
  Uncomment the print statements to see verbose outputs on why certain primers fail the search.
  '''
  for j in range(len(sequences)):
    for k in range(len(rprimer_list)):
      # print("Verifying:",organism_list[j], "(sequence)", organism_list[k], "(rprimer)")
      # verify there is indeed a product
      if j == k and PredictPCRProduct(sequences[j], fprimer, rprimer_list[k], melting_point_model) == None: 
        # print("ERROR:", organism_list[j], "PCR has no product")
        return False
      elif j != k: 
        # verify that the rprimers only bind to specified sequence
        if not PredictPCRProduct(sequences[j], fprimer, rprimer_list[k], melting_point_model) == None:
          # print("ERROR: rprimer", organism_list[j], "binds to sequence", organism_list[k])
          return False
      # print("Verified.")
  # print("Verification complete: generated rprimers only bind to specified sequence.")
  return True

def GenerateNovelPrimers(sequences, melting_point_model):
    # print(1/0)
    # Modify this function for Task 3
    # all primer sequences should be passed out in 5' -> 3' direction
    # for each sequence add a pair of primers in a tuple
    
    # for finding the same one that exists in all:
    # find match in chicken and beef
    pork_seq = sequences[0]
    chicken_seq = sequences[1]
    beef_seq = sequences[2]

    print("Searching for fprimers")
    i = 0  
    foo = 18
    while i < len(chicken_seq) - foo:
      i += 1
      if i % 100 == 0: print("Searching:", i)
      cb_score, cb_loc, A = alignment.local_align(chicken_seq[i:i+foo], beef_seq)
      if cb_score < .9 * foo * 10:
        continue
      cp_score, cp_loc, A = alignment.local_align(chicken_seq[i:i+foo], pork_seq)
      if cp_score < .9 * foo * 10:
        continue
      forward_features = CalculatePrimerFeatures(ReverseCompliment(chicken_seq[i:i+foo]))
      fT = melting_point_model.predict([forward_features])[0]
      if abs(fT - 60) > 3:
        continue
      break
    if i == len(chicken_seq) - foo: 
        print("Fprimer not found.")
        return
    else:
        print("Fprimer found! Scores:", cb_score, cp_score)

    # now we got our working fprimer

    chicken_start_i = i
    beef_start_i = cb_loc[1] - foo
    pork_start_i = cp_loc[1] - foo
    
    beef_end_i = beef_start_i + foo
    chicken_end_i = chicken_start_i + foo
    pork_end_i = pork_start_i + foo
  
    fprimer = (chicken_seq[chicken_start_i: chicken_end_i])
    
    # now we want Chicken's DNA to be: 50, Pork: 150, Beef: 250
    
    #select a location 250 away from optloc for the reverse primer
    pork_rprimer_start = pork_end_i + 50
    pork_rprimer = ReverseCompliment(pork_seq[pork_rprimer_start:pork_rprimer_start+foo])
    
    beef_rprimer_start = beef_end_i + 150
    beef_rprimer = ReverseCompliment(beef_seq[beef_rprimer_start:beef_rprimer_start+foo])
    
    chicken_rprimer_start = chicken_end_i + 250
    chicken_rprimer = ReverseCompliment(chicken_seq[chicken_rprimer_start:chicken_rprimer_start+foo])
    
    rprimer_list = [pork_rprimer, chicken_rprimer, beef_rprimer]
    organism_list = ["pork", "chicken", "beef"]
    
    print("Searching for rprimers")
    count = 0
    error_exists = True
    while error_exists:
      count += 1
      if pork_rprimer_start > len(pork_seq): # avoid index out of bound
        print("Out of bound right away. Found nothing.")
        break
      if count % 10 == 0: print("Searching:", count)
      if verify_PCR_exists(sequences, rprimer_list, organism_list, fprimer, melting_point_model):
        error_exists = False
      else:
        if PredictPCRProduct(pork_seq, fprimer, pork_rprimer, melting_point_model) == None:
          pork_rprimer_start += 1
          pork_rprimer = ReverseCompliment(pork_seq[pork_rprimer_start:pork_rprimer_start+foo])
        if PredictPCRProduct(beef_seq, fprimer, beef_rprimer, melting_point_model) == None:
          beef_rprimer_start += 1
          beef_rprimer = ReverseCompliment(beef_seq[beef_rprimer_start:beef_rprimer_start+foo])
        if PredictPCRProduct(chicken_seq, fprimer, chicken_rprimer, melting_point_model) == None:
          chicken_rprimer_start += 1
          chicken_rprimer = ReverseCompliment(chicken_seq[chicken_rprimer_start:chicken_rprimer_start+foo])
        rprimer_list = [pork_rprimer, chicken_rprimer, beef_rprimer]
    # IF THIS WORKS MAKE SURE LENGTHS ARE DIFFERENT!
    print("Rprimer found! Indices of rprimers (p,c,b):", pork_rprimer_start, chicken_rprimer_start, beef_rprimer_start)

    return [[fprimer, pork_rprimer], [fprimer, chicken_rprimer], [fprimer, beef_rprimer]]