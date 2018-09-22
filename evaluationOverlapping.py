# -*- coding: utf-8 -*-

"""	

Evaluate the motif finder: overlapping sites and running time. 
    python evaluationOverlapping.py
    @Qingqing Zhang, Zheyi Zhu
    @Chuankai Zhao, czhao37@illinois.edu

"""

import numpy as np

def getSiteFile(filename):
	f = open(filename)
	sites = []
	for line in f:
		line = line[:-1]
		sites.append(line)
	return sites

def getPredictSites(filename):
	f = open(filename)
	predictSites = []
	for line in f:
		if not line.startswith(">"):
			line = line[:-1]
			predictSites.append(line)
	return predictSites

def countOverlap(f1, icpc, ml, sl, sc, i, num):
	filePath = "datasets/" + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i) + "/"
	siteFile = filePath + "sites.txt"
	sites = getSiteFile(siteFile)
	length = len(sites)
	count = 0
	for j in range(num):
		predictFile = filePath + "predictedsites" + "_" + str(j) + ".txt"
		predictSites = getPredictSites(predictFile)
		for p in range(length):
			if sites[p] == predictSites[p]:
				count += 1
	count = count*1.0 / num
	return count 


def getRunningTime(filename):
	f = open(filename, "r+")
	t = f.read().rstrip('\n')
	h, m, s = [float(i) for i in t.split(":")]
	return 3600*h + 60*m + s
	f.close()

def calculateRunTime(f2, icpc, ml, sl, sc, i, num):
	time = 0
	filePath = "datasets/" + "dataset_" + str(icpc) + "_" + str(ml) + "_" + str(sl) + "_" + str(sc) + "_" + str(i) + "/"
	for j in range(num):
		runTimeFile = filePath + "running_time" + "_" + str(j) + ".txt"
		time += getRunningTime(runTimeFile)
	time = time / num	
	return time	


def main():
	icpcList = [1, 1.5, 2]
	mlList = [6, 7, 8]
	scList = [5, 10, 20]
	ICPC = 2
	ML = 8
	SL = 500
	SC = 10
	num = 10

	f1 = open("averageOverlap.txt", "wb+")
	f2 = open("averageRunningTime.txt", "wb+")
	for icpc in icpcList:
		total_count = []
		total_time = []
		c = 0
		t = 0
		for i in range(num):
			c = countOverlap(f1, icpc, ML, SL, SC, i, num)
			t = calculateRunTime(f2, icpc, ML, SL, SC, i, num)
			total_count.append(c)
			total_time.append(t) 
		countMean = np.mean(total_count)
		countSTD = np.std(total_count)
		timeMean = np.mean(total_time)
		timeSTD = np.std(total_time)

		f1.write("ICPC = " + str(icpc) + ", ML = " + str(ML) + ", SL = " + str(SL) + ", SC = " + str(SC) + ", Overlap Count = " + str(countMean) + ", STD = " + str(countSTD) + "\n")
		f2.write("ICPC = " + str(icpc) + ", ML = " + str(ML) + ", SL = " + str(SL) + ", SC = " + str(SC) + ", Running Time = " + str(timeMean) + ", STD = " + str(timeSTD) + "\n")
	
	for ml in mlList:
		total_count = []
		total_time = []
		c = 0
		t = 0
		for i in range(num):
			c = countOverlap(f1, ICPC, ml, SL, SC, i, num)
			t = calculateRunTime(f2, ICPC, ml, SL, SC, i, num)
			total_count.append(c)
			total_time.append(t) 
		countMean = np.mean(total_count)
		countSTD = np.std(total_count)
		timeMean = np.mean(total_time)
		timeSTD = np.std(total_time)

		f1.write("ICPC = " + str(ICPC) + ", ML = " + str(ml) + ", SL = " + str(SL) + ", SC = " + str(SC) + ", Overlap Count = " + str(countMean) + ", STD = " + str(countSTD) + "\n")
		f2.write("ICPC = " + str(ICPC) + ", ML = " + str(ml) + ", SL = " + str(SL) + ", SC = " + str(SC) + ", Running Time = " + str(timeMean) + ", STD = " + str(timeSTD) +  "\n")

	for sc in scList:
		total_count = []
		total_time = []
		c = 0
		t = 0
		for i in range(num):
			c = countOverlap(f1, ICPC, ML, SL, sc, i, num)
			t = calculateRunTime(f2, ICPC, ML, SL, sc, i, num)
			total_count.append(c)
			total_time.append(t) 
		countMean = np.mean(total_count)
		countSTD = np.std(total_count)
		timeMean = np.mean(total_time)
		timeSTD = np.std(total_time)

		f1.write("ICPC = " + str(ICPC) + ", ML = " + str(ML) + ", SL = " + str(SL) + ", SC = " + str(sc) + ", Overlap Count = " + str(countMean) + ", STD = " + str(countSTD) + "\n")
		f2.write("ICPC = " + str(ICPC) + ", ML = " + str(ML) + ", SL = " + str(SL) + ", SC = " + str(sc) + ", Running Time = " + str(timeMean) + ", STD = " + str(timeSTD) + "\n")

	f1.close()
	f2.close()


if __name__ == '__main__':
	main()
