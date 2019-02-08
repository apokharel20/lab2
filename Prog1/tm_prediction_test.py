# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
import time
import alignment
import copy
# Your task is to *accurately* predict the primer melting points using machine 
# learning based on the sequence of the primer.

# Load the primers and their melting points.

def LoadFastA(path):
    infile = open(path, 'r')
    seq = ""
    infile.readline()
    for line in infile:
        seq += line[:-1]
    return seq
    
if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
   from sklearn.metrics import r2_score
   from sklearn.ensemble import RandomForestRegressor
   import importlib
   import os
   import sys
   test_lib_files = os.listdir("test_modules")
   
   
   print("Running Task 1:")
   
   sys.path.append("test_modules")
   
   
   results = "module_name\tCrossvalidation_R2\tTesting_Set_R2\tFeatureCalcTime\tFeatureCount\n"
   for test_lib_file in test_lib_files:
       if "pycache" not in test_lib_file:    
           test_module = importlib.import_module(test_lib_file[:-3], package=None)
           print ("Loading module", test_lib_file)
           infile = open("training_primers.txt", 'r')
           infile.readline() # don't load headers
           primers = []
           melting_points = []
           features = []
           st = time.time()
           for line in infile:
               Line = line.split()
               primers.append(Line[0])
               melting_points.append(float(Line[1]))
               # calculate features
               features.append(test_module.CalculatePrimerFeatures(Line[0]))
           feat_time = (time.time()-st)/(len(features)/1000)
           # cross validation
           how_many_folds = 10 
           predictions = []
           truth = []
           my_len = len(features[-1])
            
           for fold in range(how_many_folds):
                #print ("Calculating Fold",fold)
                training_features = []
                training_outcomes = []
                testing_features = []
                testing_outcomes = []
                for c in range(len(melting_points)):
                    if c % how_many_folds == fold:
                        # put this one in testing data
                        testing_features.append(features[c])
                        testing_outcomes.append(melting_points[c])
                    else:
                        # put this one in training data
                        training_features.append(features[c])
                        training_outcomes.append(melting_points[c])
                # train the model
                
                rf = RandomForestRegressor(n_estimators = 200)
                rf.fit(training_features, training_outcomes)
                fold_predictions = rf.predict(testing_features)
                truth += testing_outcomes
                predictions += list(fold_predictions)
           
           
           
           if os.path.exists("testing_primers.txt"):
               rf = RandomForestRegressor(n_estimators = 200)
               rf.fit(features, melting_points)
               infile = open("testing_primers.txt", 'r')
               infile.readline() # don't load headers
               primers = []
               melting_points = []
               features = []
               for line in infile:
                   Line = line.split()
                   primers.append(Line[0])
                   melting_points.append(float(Line[1]))
                   # calculate features
                   features.append(test_module.CalculatePrimerFeatures(Line[0]))
               final_preds = rf.predict(features)
               
               results += "%s\t%f\t%f\t%f\t%d\n" % (test_lib_file[:-3], r2_score(truth, predictions), r2_score(melting_points, final_preds), feat_time, my_len)
           else:
               results += "%s\t%f\tN/A\t%f\t%d\n" % (test_lib_file[:-3], r2_score(truth, predictions), feat_time,my_len)
   print("Task 1 Results:\n")
   print(results)



   print()
   print("Task 2 Testing:")
   infile= open("PCR_product_test_cases.txt", 'r')
   infile.readline()
   test_pcr_conditions = []
   for line in infile:
       Line = line.split()
       #print(Line)
       test_pcr_conditions.append(copy.deepcopy(Line))
       
   for test_lib_file in test_lib_files:
       if "pycache" not in test_lib_file:    
           test_module = importlib.import_module(test_lib_file[:-3], package=None)
           print ("Loading module", test_lib_file)
           infile = open("training_primers.txt", 'r')
           infile.readline() # don't load headers
           primers = []
           melting_points = []
           features = []
           st = time.time()
           for line in infile:
               Line = line.split()
               primers.append(Line[0])
               melting_points.append(float(Line[1]))
               # calculate features
               features.append(test_module.CalculatePrimerFeatures(Line[0]))
           rf = RandomForestRegressor(n_estimators = 200)
           rf.fit(features, melting_points)
           # generate new sequences for testing
           c = 0
           
           for t in test_pcr_conditions:
               true_product = t[3]
               primer1 = t[1]
               primer2 = t[2]
               full_seq = t[0]
               print("Case:", c,"\tTruth:", true_product,"\tPrediction:", test_module.PredictPCRProduct(full_seq,primer1,primer2, rf))
               c += 1
   print("\n")
   print("Task 3 Run:")
   pork_cytob = LoadFastA('sus_cytochrome_b.fa')
   chicken_cytob = LoadFastA('gallus_cytochrome_b.fa')
   beef_cytob = LoadFastA('bos_taurus_cytochrome_b.fa')
   
   pork_mito = LoadFastA('sus_scrofa_domestica_mtDNA.fa')
   chicken_mito = LoadFastA('gallus_gallus_mtDNA.fasta')
   beef_mito = LoadFastA('bos_taurus_mtDNA.fasta')
   
   for test_lib_file in test_lib_files:
       if "pycache" not in test_lib_file:    
           test_module = importlib.import_module(test_lib_file[:-3], package=None)
           print ("Loading module", test_lib_file)
           infile = open("training_primers.txt", 'r')
           infile.readline() # don't load headers
           primers = []
           melting_points = []
           features = []
           st = time.time()
           for line in infile:
               Line = line.split()
               primers.append(Line[0])
               melting_points.append(float(Line[1]))
               # calculate features
               features.append(test_module.CalculatePrimerFeatures(Line[0]))
           rf = RandomForestRegressor(n_estimators = 200)
           rf.fit(features, melting_points)
           # generate new sequences for testing
           
           sequences = [pork_cytob, chicken_cytob, beef_cytob]
           
           print("Designed Primers:", test_module.GenerateNovelPrimers(sequences, rf))
           

# Calculate features

