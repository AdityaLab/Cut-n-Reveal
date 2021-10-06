import os
import sys
import networkx as nx
import ast
import matplotlib
import numpy as np
matplotlib.use('AGG')
import matplotlib.pyplot as plt
def plot_result(data_dir, hurricane, k, save_dir):
    """
    plot the network behind the time series for each cut point, highlight the important time series in the explanation accordingly.
    also plot the time series data and the segmentation
    """
    k = int(k)
    save_dir += hurricane + '/' + str(k) +'/'
    #with open(data_dir + hurricane + '.outages_data') as f:
    #    line = f.readline()
    #counties = line.strip().split('\t')
    #with open(data_dir + 'TICC/' + hurricane + '/interp_data.csv') as f:
    #    lines = f.readlines()
    with open(data_dir + hurricane + '_interp.csv')as f:
    #with open(data_dir + 'hurricane_matthew_60min_sampling_interval.csv') as f:
    #with open(data_dir + 'Matthew.csv') as f:
        lines = f.readlines()
    counties = lines[0].strip().split(',')
    lines = lines[1:]
    data = []
    for line in lines:
        data.append([float(x) for x in line.strip().split(',')])
    data = np.array(data)
    l = len(lines)
    print 'reading the network'
    G = nx.read_adjlist(data_dir + hurricane + '.adjlist', delimiter = '\t')
    print 'reading the segmentation result'
    with open(save_dir + 'S_' + str(k) + '.txt') as f:
        line = f.readline()
    S = ast.literal_eval(line.strip())
    print 'reading the explanation'
    E = np.load(save_dir + 'E_' + str(k) + '.txt.npy')
    #E = (E - E.min(0)) / (E.max(0) - E.min(0))

    print 'plotting the segmentation'
    for i in range(len(counties)):
        plt.plot(range(l), data[:, i])
    for s in S:
        plt.axvline(x=float(s))
    plt.savefig(save_dir + '/S_' + str(k) + '.pdf')

    plt.clf()
    print 'plotting the graphs'
    node_labels = {}
    nodes = G.nodes()
    #for i in range(k):
    for i in range(len(S)):
        node_size = []
        for n in nodes:
            index = counties.index(n)
            w = E[index, i] * 300
            node_size.append(w)
        nx.draw_networkx(G, nodelist = nodes, node_size = node_size, font_size = 5)
        plt.savefig(save_dir + '/Egraph_' + str(k) + '_cut' + str(i) + '.pdf')
        plt.clf()
if __name__ == '__main__':
    data_dir = './data/'
    #hurricane = 'Harvey'
    hurricane = 'syn2'
    k = sys.argv[1]
    #save_dir = '../result/'
    #save_dir = './HurricaneTimeseriesSegmentation/source/baseline/SubspaceClustering/result/'
    save_dir = './result/'
    plot_result(data_dir, hurricane, k, save_dir)
