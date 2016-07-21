# -*- coding: utf-8 -*-
"""
Created on Mon May 30 17:36:01 2016

@author: colinrathbun
"""

import math
import numpy as np
import pandas as pd
import itertools as it
import csv


def calc_orthog_index(data, cmpd1, cmpd2, mut1, mut2):
    value = math.log10(abs(float(data[mut1, cmpd1]) / float(data[mut1, cmpd2]))) * math.log10(abs(float(data[mut2, cmpd2]) / float(data[mut2, cmpd1])))
    if np.isnan(value):
        value = -1
    return value


def import_raw_data(path):
    rawData = list(csv.reader(open(path, "r"), delimiter=','))
    fullData = np.genfromtxt(path, delimiter=',')
    #fullData = np.array(rawData, dtype=np.float)
    data = np.delete(fullData, 0, 0)
    data = np.delete(data, 0, 1)
    #Remove outliers and set all values below 1E3 to 1E3
    for flux_value in np.nditer(data, op_flags=['readwrite']):
        if flux_value < 1000:
            flux_value[...] = 1000 #using the elipsis will actually set the value in the array
    compounds = rawData[0][1:]
    mutants = []
    for line in rawData:
        mutants.append(line[0])
    mutants = mutants[1:]
    return [data, compounds, mutants]

def compute_orthog_and_rank(data):
    compound_combinations = it.combinations(range(data.shape[1]), 2)
    mutant_combinations = it.combinations(range(data.shape[0]), 2)
    all_combinations = it.product(compound_combinations, mutant_combinations)
    
    num_combinations = math.ceil(data.shape[1] * data.shape[0] * (data.shape[1] - 1) * (data.shape[0] - 1) / 4)
    #orthog_results = np.zeros((num_combinations, 5), dtype=[('orthog', float), ('cmpd1', int), ('cmpd2', int), ('mut1', int), ('mut2', int)])
    orthog_results = np.zeros((num_combinations, 5))
    
    c = 0
    for combination in all_combinations:
        orthog = calc_orthog_index(data, combination[0][0], combination[0][1], 
                                   combination[1][0], combination[1][1])
        orthog_results[c, 0] = orthog
        orthog_results[c, 1] = combination[0][0]
        orthog_results[c, 2] = combination[0][1]
        orthog_results[c, 3] = combination[1][0]
        orthog_results[c, 4] = combination[1][1]
        c += 1
    
    orthog_results.view('f8,f8,f8,f8,f8').sort(order=['f0'], axis=0)
    orthog_results = np.flipud(orthog_results)
    return orthog_results

#orthog_results = np.sort(orthog_results, order='f0')
#orthog_results = sorted(orthog_results, key=lambda i: i[0])
#orthog_results.reverse()

def get_readable_list(rawData, results, compounds, mutants, number_of_top):
    c = 0    
    working_list = []
    while c < number_of_top:
        c1 = int(results[c, 1])
        c2 = int(results[c, 2])
        m1 = int(results[c, 3])
        m2 = int(results[c, 4])
        if rawData[m1, c1] > rawData[m1, c2]:
            working_list.append([results[c][0], compounds[c1], mutants[m1], rawData[m1,c1], rawData[m1,c1]/rawData[m1,c2], compounds[c2], mutants[m2], rawData[m2,c2], rawData[m2,c2]/rawData[m2,c1]])
        else:
            working_list.append([results[c][0], compounds[c1], mutants[m2], rawData[m2,c1], rawData[m2,c1]/rawData[m2,c2], compounds[c2], mutants[m1], rawData[m1,c2], rawData[m1,c2]/rawData[m1,c1]])
        c += 1
    return pd.DataFrame(np.array(working_list), index=list(range(number_of_top)), columns=['orthog', 'cmpd1', 'mut1', 'flux1', 'reso1', 'cmpd2', 'mut2', 'flux2', 'reso2'])

def readable_list_from_path(path, number_of_top):
    rawDatas = import_raw_data(path)
    data = rawDatas[0]
    compounds = rawDatas[1]
    mutants = rawDatas[2]
    unreadable = compute_orthog_and_rank(data)
    return get_readable_list(data, unreadable, compounds, mutants, number_of_top)
    
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


    
        