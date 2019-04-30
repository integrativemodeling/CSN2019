import sys, os
import math
import numpy as np


def get_scores_distribution(scores, nbins, scorecolumn, hist_ofile):
    
    H, xedge = np.histogram(scores, bins=nbins)
    f1=open(hist_ofile, 'w+')
    
    for i in range(nbins):
        print >>f1, xedge[i], H[i]
    return

scores = []
with open(sys.argv[1], 'r') as lines:
    for line in lines:
        line = line.strip().split(" ")[0]
        line = float(line)
        scores.append(float(line))

get_scores_distribution(scores, 50, 0, sys.argv[2])
