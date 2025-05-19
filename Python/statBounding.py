'''
Sample from normal distribution, and calculate bounding leakage.
'''
import datetime
import itertools
import math
import numpy
import random
import scipy.special

SAMPLES = (2,3,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000)
KAPPA = (1,2,3,4,5,6)
COUNT = 10000


ssStat = {}

def calcStat():
    filePath = "./Python/Output/StatBounding.txt"
    print(f'Start output to {filePath}')
    with open(filePath, 'w') as f:
        f.write('Samples\tKappa\tCount\tMean\tDeviation\n')
        for k in KAPPA:
            f.write(f'0\t{k}\t{COUNT}\t{1 - scipy.special.erf(k/math.sqrt(2))}\t0\n')
        f.flush()
        for samples in SAMPLES:
            print(f'Start samples={samples} at {datetime.datetime.now()}')
            sStat = {k:[0]*3 for k in KAPPA}
            ssStat[samples] = sStat
            for i in range(COUNT):
                sSample = [random.gauss() for j in range(samples)]
                mu = numpy.average(sSample)
                sigma = numpy.std(sSample)
                for k in KAPPA:
                    err = 1 - (scipy.special.erf((k + mu)/sigma/math.sqrt(2)) + scipy.special.erf(abs((k - mu)/sigma/math.sqrt(2))))/2
                    stat = sStat[k]
                    stat[0] += 1
                    stat[1] += err
                    stat[2] += err**2
            for k in KAPPA:
                stat = sStat[k]
                stat[1] /= stat[0]
                f.write(f'{samples}\t{k}\t{stat[0]}\t{stat[1]}\t{math.sqrt(stat[2]/stat[0] - stat[1]**2)}\n')
            f.flush()



if __name__ == '__main__':
    calcStat()