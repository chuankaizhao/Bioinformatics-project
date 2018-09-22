"""

Calculate relative entropy between motif.txt and predictedmotif.txt and evaluate the motif finder.
    python relativeEntropy.py

@Chuankai Zhao, czhao37@illinois.edu
@Zheyi Zhu, Qingqing Zhang

"""

import os
import math
import sys
import random
import numpy as np

def get_profile(dataset_dir, filename):

    f1 = open(dataset_dir + filename, "r")
    line = f1.readlines()
    line = line[1:-1]

    motif = []
    for i in range(len(line)):
        if ">" or "<" not in line[i]:
            motif.append(line[i].rstrip().split('\t'))

    for i in range(len(motif)):
        for j in range(len(motif[i])):
            motif[i][j] = int(motif[i][j])
 
    return motif

# generate the position weight matrix from the profile matrix
def get_pwm(motif):

    pwm = []

    for i in range(len(motif)):
        pwm_column = []
        sum = np.sum(motif[i]) + 1

        # calculate the pseudocount probability (To avoid p = 0)
        for j in range(len(motif[0])):
            weight = (motif[i][j] + 0.25)/sum
            pwm_column.append(weight)
        pwm.append(pwm_column)

    return pwm

# get information content
def get_info_content(pwm):

    W = 0.

    for i in range(len(pwm)):
        for j in range(len(pwm[0])):
            W = W + pwm[i][j] * np.log2( pwm[i][j] * 4.0 )

    return W

# calculate the information content from the position weight matrix
def get_relative_entropy(pwm_m, pwm_p):

    RE = 0.
    
    for i in range(len(pwm_p)):
        for j in range(len(pwm_p[0])):
            RE = RE + pwm_p[i][j] * np.log2( pwm_p[i][j] / pwm_m[i][j] )
    
    return RE

# write the predicted motif, predicted sites, running time and running information into files
def write(icpc, ml, sl, sc, re, re_d):

    f1 = open('averageRelativeEntropy.txt', 'a')
    f1.write("ICPC = " + str(icpc) + ", ML = " + str(ml) + ", SL = " + str(sl) + ", SC = " + str(sc) + ", Relative Entropy = " + str(re) + ", Standard Error = " + str(re_d) + "\n")
    f1.close()

def writeBest(icpc, ml, sl, sc, re, re_d):

    f1 = open('averageBestRelativeEntropy.txt', 'a')
    f1.write("ICPC = " + str(icpc) + ", ML = " + str(ml) + ", SL = " + str(sl) + ", SC = " + str(sc) + ", Relative Entropy = " + str(re) + ", Standard Error = " + str(re_d) + "\n")
    f1.close()

def runBest(icpc, ml, sl, sc, num_dataset, num_runs):

    datasets_directory = "./datasets/"

    res = []

    for i in range(num_dataset):

        dataset_dir = datasets_directory + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i) + "/"

        motif = get_profile(dataset_dir, "motif.txt")
        pwm_m = get_pwm(motif)

        info_cont = []
        for j in range(num_runs):
            predicted_motif = get_profile(dataset_dir, "predictedmotif" + "_" + str(j) + ".txt")
            pwm_p  = get_pwm(predicted_motif)
            info   = get_info_content(pwm_p)
            info_cont.append(info)

        num = np.argmax(info_cont)
        predicted_motif = get_profile(dataset_dir, "predictedmotif" + "_" + str(num) + ".txt")
        pwm_p  = get_pwm(predicted_motif)
        re     = get_relative_entropy(pwm_m, pwm_p)
        res.append(re)

    re = np.mean(res)
    re_d = np.std(res)
    return re, re_d


def run(icpc, ml, sl, sc, num_dataset, num_runs):
    
    datasets_directory = "./datasets/"

    res = []

    for i in range(num_dataset):
        
        dataset_dir = datasets_directory + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i) + "/"

        motif = get_profile(dataset_dir, "motif.txt")
        pwm_m = get_pwm(motif)

        for j in range(num_runs):
            predicted_motif = get_profile(dataset_dir, "predictedmotif" + "_" + str(j) + ".txt")
            pwm_p  = get_pwm(predicted_motif)
            re     = get_relative_entropy(pwm_m, pwm_p)     
            res.append(re) 

    re = np.mean(res)
    re_d = np.std(res)
    return re, re_d  

if __name__ == '__main__':

    ICPC = 2
    ML = 8
    SL = 500
    SC = 10
    ICPC_list = [1,1.5,2]
    ML_list = [6,7,8]
    SC_list = [5,10,20]
    num_dataset = 10
    num_runs    = 10

    for icpc in ICPC_list:
        ml = ML
        sl = SL
        sc = SC
        re, re_d = run(icpc, ml, sl, sc, num_dataset, num_runs)
        write(icpc, ml, sl, sc, re, re_d)
        re, re_d = runBest(icpc, ml, sl, sc, num_dataset, num_runs)
        writeBest(icpc, ml, sl, sc, re, re_d)

    for ml in ML_list:
        icpc = ICPC
        sl   = SL
        sc   = SC
        re, re_d = run(icpc, ml, sl, sc, num_dataset, num_runs)
        write(icpc, ml, sl, sc, re, re_d)
        re, re_d = runBest(icpc, ml, sl, sc, num_dataset, num_runs)
        writeBest(icpc, ml, sl, sc, re, re_d)

    for sc in SC_list:
        icpc = ICPC
        sl   = SL
        ml   = ML
        re, re_d = run(icpc, ml, sl, sc, num_dataset, num_runs)
        write(icpc, ml, sl, sc, re, re_d)
        re, re_d = runBest(icpc, ml, sl, sc, num_dataset, num_runs)
        writeBest(icpc, ml, sl, sc, re, re_d)
