# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 11:15:31 2016

@author: colinrathbun
"""

import math
import numpy as np
import pandas as pd
import itertools as it
import csv
#from progress.bar import Bar
    
def sets_are_nodes3(setlist):
    edgeList = []
    nodeList = []
    for member in setlist:
        edgeList.append(tuple(sorted([member[1],member[2],member[4],member[5]])))
        nodeList.append((member[1],member[2]))
        nodeList.append((member[4],member[5]))
    result1 = result2 = False
    if len(set(edgeList)) == len(setlist):
        result1 = True
    if len(set(nodeList)) == len(setlist):
        result2 = True
    return result1 and result2

def get_triangle_list(pdset):
    
    #setup a progress bar
    #bar = Bar('Iterating', max=len(pdset))
    denom = len(pdset)
    numer = 0    
    
    npSet = pdset.values
    triangles = []
    for set1 in npSet:
        for set2 in npSet:
            for set3 in npSet:
                setlist = (set1,set2,set3)
#                setlist = sorted(((set1[1],set1[2],set1[4],set1[5]),(set2[1],set2[2],set2[4],set2[5]),(set3[1],set3[2],set3[4],set3[5]))
                if sets_are_nodes3(setlist):
                    triangles.append(setlist)
        numer += 1
        print("Progress: "+str(numer / denom*100)+"% Found "+str(len(triangles)))
    return triangles

def write_out_triangles(trianglesList, path):
	with open(path, 'wb') as f:
		np.savetxt(f,np.reshape(np.array(trianglesList), (-1,7)), fmt='%s, %s, %s, %s, %s, %s, %s')

