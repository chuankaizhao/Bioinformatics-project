"""

Running motif finding using gibbs sampling algorithm.
    python motifFindingGibbs.py icpc ML SL SC datasetID N_runs 
    eg. python motifFindingGibbs.py 2 8 500 10 0 10

@Chuankai Zhao, czhao37@illinois.edu
@Zheyi Zhu, Qingqing Zhang

"""

import os
import math
import sys
import random
import numpy as np
from datetime import datetime

# generate random motif as the initial motif finding solution
def random_motif(seqs, ML):

    pos = []

    for i in range(len(seqs)):
        pos.append(random.randrange(0,len(seqs[0])-ML))

    aligned_patterns = get_aligned_patterns(seqs, pos, ML)
    motif = get_profile(aligned_patterns)

    return pos, motif

# generate the aligned motif patterns 
def get_aligned_patterns(seqs, pos, ML):

    aligned_patterns = []

    for i in range(len(seqs)):
        aligned_patterns.append(seqs[i][pos[i]:pos[i]+ML])

    return aligned_patterns

# generate the profile matrix from the aligned motif patterns
def get_profile(aligned_patterns):

    motif = []

    for i in range(len(aligned_patterns[0])):
        counts = np.zeros((4),dtype=int)

        for j in range(len(aligned_patterns)):
            if aligned_patterns[j][i] == 'A': counts[0] = counts[0] + 1
            if aligned_patterns[j][i] == 'C': counts[1] = counts[1] + 1
            if aligned_patterns[j][i] == 'G': counts[2] = counts[2] + 1
            if aligned_patterns[j][i] == 'T': counts[3] = counts[3] + 1

        for k in range(4):
            counts[k] = int(counts[k])

        motif.append(counts)

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

# calculate the information content from the position weight matrix
def get_info_content(pwm):

    W = 0.

    for i in range(len(pwm)):
        for j in range(len(pwm[0])):
            W = W + pwm[i][j] * np.log2( pwm[i][j] * 4.0 )

    return W

# implement the gibbs sampling algorithm to do the motif finding. 
def motif_finding(seqs, ML):

    # generate random motif as initial solution
    pos, motif = random_motif(seqs, ML)

    # record all the samples collected from many iterations
    pos_database  = []
    info_database = []
    running_info  = []

    ITER_CONTINUE_FLAG = 1
    iter          = 0
    equal_maxinfo_N = 0

    while ITER_CONTINUE_FLAG: 
        for i in range(len(seqs)):

            # delete sequence i and generate new pwm
            seqs_rd = seqs[0:i] + seqs[i+1:]
            pos_rd  = pos[0:i]  + pos[i+1:]
            aligned_patterns = get_aligned_patterns(seqs_rd, pos_rd, ML)
            motif   = get_profile(aligned_patterns)
            pwm     = get_pwm(motif)

            # find the optimal new pos x in seq i
            
            ## calculate the probability of x from the new pwm
            prob_qx = []
            dict    = { 'A':0, 'C':1, 'G':2, 'T':3 }

            for j in range(len(seqs[0])-ML):
                prob = 1.0
                for k in range(ML):
                    prob = prob * pwm[k][dict[seqs[i][j+k]]]
                prob_qx.append(prob)
            prob_qx_norm  = prob_qx / np.sum(prob_qx)

            ## randomly generate the new pos x based on their probability distribution
            sum_prob_qx   = []
            sum      = 0.0

            for q in range(len(prob_qx_norm)):
                sum  =  sum + prob_qx_norm[q]
                sum_prob_qx.append(sum)

            rand   = random.random()

            for m in range(len(prob_qx_norm)):
                if m == 0:
                    if rand < sum_prob_qx[m]: max_pos = m
                if m > 0:
                    if sum_prob_qx[m-1] < rand < sum_prob_qx[m]: max_pos = m
                if m == len(prob_qx_norm) - 1:
                    if sum_prob_qx[m] < rand: max_pos = m
            pos = pos[0:i] + [max_pos] + pos[i+1:]

        # generate the new motif after new pos x in seq i is found
        aligned_patterns = get_aligned_patterns(seqs, pos, ML)
        motif = get_profile(aligned_patterns)
        pwm     = get_pwm(motif)

        # calculate the infomation content for this iteration
        info    = get_info_content(pwm)
        pos_database.append(pos)
        info_database.append(info) 

        # judge whether to continue iterations or not
        if iter > 999:
            if np.max(info_database) == max_info:
                equal_maxinfo_N = equal_maxinfo_N + 1
            if np.max(info_database) != max_info: 
                equal_maxinfo_N = 0
            if equal_maxinfo_N > 499: 
                ITER_CONTINUE_FLAG = 0

        # record the max_info
        max_info = np.max(info_database)
        iter = iter + 1
        running_info.append([iter, equal_maxinfo_N, info, max_info])

    # generate the predicted motif and sites
    max_info_pos = np.argmax(info_database)
    max_pos      = pos_database[max_info_pos]
    aligned_patterns = get_aligned_patterns(seqs, max_pos, ML)
    motif    = get_profile(aligned_patterns)

    base = ['A','C','G','T']
    final_motif = ''
    for i in range(len(motif)):
        final_motif = final_motif + base[np.argmax(motif[i])]

    return pos, motif, final_motif, running_info

# write the predicted motif, predicted sites, running time and running information into files
def write_to_file(dataset_dir, pos, motif, ml_ri, eclapse, final_motif, runiter, running_info):
    
    sites_dir = dataset_dir + '/predictedsites_' + str(runiter) + '.txt'
    f = open(sites_dir, 'w')
    for i in range(len(pos)):
        f.write('>site' + str(i+1) + '\n')
        f.write(str(pos[i]))
        f.write("\n")
    f.close()
    
    motif_dir = dataset_dir + '/predictedmotif_' + str(runiter) + '.txt'
    f = open(motif_dir, 'w')
    f.write('>motif\t' + final_motif + '\t' + str(ml_ri) + '\n')
    for i in range(len(motif)):
        f.write('\t'.join(map(str, motif[i])))
        f.write("\n")
    f.write("<\n")
    f.close()

    runtime_dir = dataset_dir + '/running_time_' + str(runiter) + '.txt'
    f = open(runtime_dir, 'w')
    f.write(str(eclapse) + '\n')
    f.close()
  
    running_info_dir = dataset_dir + '/running_info_' + str(runiter) + '.txt'
    f = open(running_info_dir, 'w')
    for i in range(len(running_info)):
        f.write('\t'.join(map(str, running_info[i])))
        f.write('\n')
    f.close() 

# read the sequences and motif length, and run the motif finding process
def run(icpc, ml, sl, sc, id, num):

    datasets_directory = "./datasets/"
    dataset_dir = datasets_directory + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(id)

    ## read in sequences
    f1 = open(dataset_dir+'/sequences.fa', 'r')
    line = f1.readlines()
    seqs = []
    for i in range(len(line)):
        if '>' not in line[i]:
            seqs.append(line[i].replace("\n",''))

    ## read in ML
    f2 = open(dataset_dir+'/motiflength.txt', 'r')
    ML_ri = int(f2.readlines()[0])

    ## motif finding & test time
    for runiter in range(num):
        time_begin = datetime.now()
        pos, motif, final_motif, running_info = motif_finding(seqs, ML_ri)
        time_finish = datetime.now()
        eclapse = time_finish - time_begin
        ## write results to file
        write_to_file(dataset_dir, pos, motif, ML_ri, eclapse, final_motif, runiter, running_info)

if __name__ == '__main__':
    icpc = sys.argv[1]
    ml   = int(sys.argv[2])
    sl   = int(sys.argv[3])
    sc   = int(sys.argv[4])
    id   = int(sys.argv[5])
    num  = int(sys.argv[6])
    
    run(icpc, ml, sl, sc, id, num)
