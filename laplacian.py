# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import networkx as nx
import pandas as pd
import sys
from scipy.io import savemat
from scipy.io import loadmat
import numpy as np

if __name__ == '__main__':
    fileName = sys.argv[1] #'../us_state_covid.csv'
    adjlist = sys.argv[2] #'../us_states.adjlist'
    savefile= sys.argv[3] #'../us_states'
    data = pd.read_csv(fileName)
    counties = list(data.columns.values)
    #print(len(counties))
    G = nx.read_adjlist(adjlist,delimiter=',')
    print(G.edges)
    L = nx.laplacian_matrix(G,counties)
    np.savetxt(savefile+'.csv',L.toarray(),delimiter=',')
    #print(L)
    #savemat(savefile,{'L':L})
    #mat_file = loadmat(savefile)
    #print(mat_file['L'])

