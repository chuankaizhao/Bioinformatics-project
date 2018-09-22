"""

Generate the synthetic dataset used for the motif finding.
    python generateDataset.py

@Chuankai Zhao, czhao37@illinois.edu 
@Zheyi Zhu, Qingqing Zhang

"""

import numpy as np
import random
import os

# Genearte random sequences
def generate_sequences(SL, SC):
	rd_seq = []
	for i in range(SC):
		rd_seq.append(''.join(random.choice("ACGT") for i in range(SL)))
	return rd_seq

# Generate ICPC Matrix with given ICPC value
def genarate_ICPC_Matrix(ICPC):
	ICPC_Matrix = []

	if ICPC == 1.0:
		p = 0.8105
	elif ICPC == 1.5:
		p = 0.9245  
	elif ICPC == 2.0:
		p = 1.0
	q = (1-p)/3
	
	for i in range(0, 4):
		newLine = []
		for j in range(0, 4):
			if j == i:
				newLine.append(p)
			else:
				newLine.append(q)
		ICPC_Matrix.append(newLine)
	return ICPC_Matrix

# Generate random motif_matrix and motif_string 
def generate_motifs(ICPC, ML, SC):
	# @matrix_matrix: sum of ACGT values at each position with length ML
 	# @matrix_string: collection of random motif strings with size SC 

	BASE = ["A", "C", "G", "T"]
	ICPC_Matrix = genarate_ICPC_Matrix(ICPC)
	motif_matrix = [[0 for m in range(4)] for n in range(ML)]
	motif_string = ["" for i in range(SC)]
	for i in range(ML):
		idx = random.randint(0,3)
		probList = ICPC_Matrix[idx]
		for j in range(SC):
			# Random generate one char with ICPC distribution
			motif_char_list = np.random.choice(BASE, 1, p = probList) 
			motif_char = motif_char_list[0]
			motif_string[j] += motif_char

			if motif_char == "A":
				motif_matrix[i][0] += 1		
			elif motif_char == "C":
				motif_matrix[i][1] += 1	
			elif motif_char == "G":
				motif_matrix[i][2] += 1	
			elif motif_char == "T":	
				motif_matrix[i][3] += 1

	return motif_matrix, motif_string

def generate_binding_site(SC, SL, ML):
	binding_sites = []
	for i in range(SC):
		binding_sites.append(random.randint(0, SL-ML))
	return binding_sites

def generate_planted_sequences(ICPC, ML, SL, SC):
	new_rd_seq = []
	rd_seq = generate_sequences(SL, SC)
	binding_sites = generate_binding_site(SC, SL, ML)
	motif_matrix, motif_string = generate_motifs(ICPC, ML, SC)

	for i in range(SC):
		new_seq_line = ""
		start = binding_sites[i]
		new_seq_line = rd_seq[i][:start] + motif_string[i] + rd_seq[i][start+ML:]
		new_rd_seq.append(new_seq_line)
	return new_rd_seq, binding_sites, motif_matrix


def write_to_file(dataset_dir, binding_sites, new_rd_seq, motif_matrix, ML):
    
    seqs_dir = dataset_dir + '/sequences.fa'
    f = open(seqs_dir, 'w')
    for i in range(len(new_rd_seq)):
        f.write('>sequence' + str(i+1) + '\n')
        f.write(new_rd_seq[i])
        f.write("\n")
    f.close()
    
    sites_dir = dataset_dir + '/sites.txt'
    f = open(sites_dir, 'w')
    for i in range(len(binding_sites)):
        #f.write('>site' + str(i+1) + '\n')
        f.write(str(binding_sites[i]))
        f.write("\n")
    f.close()
    
    motif_dir = dataset_dir + '/motif.txt'
    f = open(motif_dir, 'w')
    f.write('>motif\t' + str(ML) + '\n')
    for i in range(len(motif_matrix)):
        f.write('\t'.join(map(str, motif_matrix[i])))
        f.write("\n")
    f.write("<\n")
    f.close()
    
    ml_dir = dataset_dir + '/motiflength.txt'
    f = open(ml_dir, 'w')
    f.write(str(ML))
    f.close() 


def main():
	ICPC = 2
	ML = 8
	SL = 500
	SC = 10
	ICPC_list = [1, 1.5]
	ML_list = [6, 7]
	SC_list = [5, 20]
	sl = SL
	num_dataset = 10

	datasets_directory = "./datasets/"
	for icpc in ICPC_list:
	    ml = ML
	    sc = SC
	    for i in range(num_dataset):            
	        dataset_directory = datasets_directory + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i)
	        if not os.path.exists(dataset_directory):                
	            os.makedirs(dataset_directory)
	            new_rd_seq, binding_sites, motif_matrix = generate_planted_sequences(icpc, ml, sl, sc)
	            write_to_file(dataset_directory, binding_sites, new_rd_seq, motif_matrix, ml)
	        
	for ml in ML_list:
	    icpc = ICPC
	    sc = SC
	    for i in range(num_dataset):            
	        dataset_directory = datasets_directory + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i)
	        if not os.path.exists(dataset_directory):                
	            os.makedirs(dataset_directory)
	            new_rd_seq, binding_sites, motif_matrix = generate_planted_sequences(icpc, ml, sl, sc)
	            write_to_file(dataset_directory, binding_sites, new_rd_seq, motif_matrix, ml)
	        
	for sc in SC_list:
	    icpc = ICPC
	    ml = ML
	    for i in range(num_dataset):            
	        dataset_directory = datasets_directory + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i)
	        if not os.path.exists(dataset_directory):                
	            os.makedirs(dataset_directory)
	            new_rd_seq, binding_sites, motif_matrix = generate_planted_sequences(icpc, ml, sl, sc)
	            write_to_file(dataset_directory, binding_sites, new_rd_seq, motif_matrix, ml)


	icpc = ICPC
	ml = ML
	sc = SC
	for i in range(num_dataset):            
	    dataset_directory = datasets_directory + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i)
	    if not os.path.exists(dataset_directory):                
	        os.makedirs(dataset_directory)
	        new_rd_seq, binding_sites, motif_matrix = generate_planted_sequences(icpc, ml, sl, sc)
	        write_to_file(dataset_directory, binding_sites, new_rd_seq, motif_matrix, ml)

if __name__ == '__main__':
	main()
