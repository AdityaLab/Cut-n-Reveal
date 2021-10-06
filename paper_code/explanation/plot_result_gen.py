import os
import sys
import networkx as nx
import ast
import random
import matplotlib
import numpy as np
matplotlib.use('AGG')
import matplotlib.pyplot as plt
def plot_result(data_dir, clique, adj_file, filename, E_file, segfile, save_dir, savefile, thres, hurricane):
    """
    plot the network behind the time series for each cut point, highlight the important time series in the explanation accordingly.
    also plot the time series data and the segmentation
    """
    with open(data_dir + filename)as f:
        lines = f.readlines()
    counties = lines[0].strip().split(',')
    lines = lines[1:]
    data = []
    for line in lines:
        data.append([float(x) for x in line.strip().split(',')])
    data = np.array(data)
    l = len(lines)
    #print 'reading the network'
    if not clique:
        G = nx.read_adjlist(data_dir + adj_file, delimiter = '\t')
    else:
        G = nx.Graph()
        for a in counties:
            for b in counties:
                if a == b:
                    continue
                G.add_edge(a, b)
    print 'reading the segmentation result'
    with open(save_dir + segfile) as f:
        line = f.readline()
    #S = ast.literal_eval(line.strip())
    S = [int(x) for x in line.strip().split(',')]
    print 'reading the explanation'
    E = np.loadtxt(save_dir + E_file)
    #E = (E - E.min(0)) / (E.max(0) - E.min(0))

    colormap = plt.cm.gist_ncar
    county_index = np.linspace(0, 0.9, len(counties))
    #random.shuffle(county_index)

    #plt.gca().set_color_cycle([colormap(i) for i in county_index])
    plt.gca().set_prop_cycle('color',[colormap(i) for i in county_index])
    print 'plotting the segmentation'
    for i in range(len(counties)):
        plt.plot(range(l), data[:, i])
    for s in S:
        plt.axvline(x=float(s), linestyle = 'dashed', color = 'k')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=15, rotation = 45)
    plt.savefig(save_dir + savefile +'.pdf')
    plt.clf()
    print 'plotting the explanations'
    segments = []
    bd = 0
    for i in range(len(S)):
        if i - 1 >= 0:
            lc = S[i - 1] + bd
        else:
            lc = 0 + bd 
        if i + 1 < len(S):
            rc = S[i + 1] - bd
        else:
            rc = len(data) - bd
        segments.append([lc, rc])
    sf = open(save_dir + savefile + '_cut' + '_impC.txt', 'wb')
    for i in range(len(segments)):
        imp_c = get_imp_c(E, i, thres) 
        sf.write('cut ' + str(i) + '\n')
        for ic in imp_c:
            sf.write(counties[ic] + '\t' + str(E[ic, i]) + '\n')
        sf.write('\n')
        seg = segments[i]
        for c in imp_c:
            #print seg[0], seg[1], c
            plt.plot(range(seg[0], seg[1]), data[seg[0]:seg[1], c])
    sf.close()
    for s in S:
        plt.axvline(x=float(s), linestyle = 'dashed', color = 'k')
    plt.savefig(save_dir + savefile +'_exp.pdf')
    plt.clf()
    for i in range(len(segments)):
        imp_c = get_imp_c(E, i, thres) 
        seg = segments[i]
        if hurricane == 'Harvey':
            ym = 120000
            xm = 300
        elif hurricane[:4]== 'Irma':
            ym = 900000
            xm = 180
        elif hurricane == 'Matthew':
            ym = 350000
            xm = 700#250
        for c in imp_c:
            #print seg[0], seg[1], c
            plt.plot(range(seg[0], seg[1]), data[seg[0]:seg[1], c], color = colormap(county_index[c]))
        for j in range(len(S)):
            s = S[j]
            if abs(j - i) <= 1: 
                plt.axvline(x=float(s), linestyle = 'dashed', color = 'k')
        plt.xlim(-10, xm)
	#plt.xlim(right=xm)
        plt.ylim(0, ym)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=15, rotation = 45)
        plt.savefig(save_dir + savefile + '_exp_' + str(i) + '.pdf')
        plt.clf()
    """
    print 'plotting the graphs'
    node_labels = {}
    nodes = G.nodes()
    for i in range(len(S)):
        node_size = []
        for n in nodes:
            index = counties.index(n)
            w = E[index, i] * 300
            node_size.append(w)
        nx.draw_networkx(G, nodelist = nodes, node_size = node_size, font_size = 5)
        plt.savefig(save_dir + savefile + '_cut' + str(i) + '.pdf')
        plt.clf()
    """
def get_imp_c(E, i, thres):
    #given the explanations, return for the ith cut point, the time series index that account for more than thres of the importance
    items = zip(range(len(E)), E[:, i])
    items = sorted(items, reverse = True, key = lambda x:x[1])
    imp_c = []
    #thres = 0.8
    s = 0.0
    for item in items:
        s += item[1]
        imp_c.append(item[0])
        if s > thres:
            break
    return imp_c
if __name__ == '__main__':
    #data_dir = '../data/non_intrusive_load_monitoring/'
    #filename = 'non_intrusive_load_monitoring_dataset_3.csv'
    #segfile = 'osc_segment_indices_lambda_1_0.5_lambda_2_1000_numiter_300.csv'
    #E_file = 'E_23.txt.npy'
    #save_dir = '../result/nilm3/'
    #non = 'normalized' 

    #data_dir = '../data/'
    #filename = 'Matthew_60min_sample.csv'
    #save_dir = '../result/Matthew_Normalized/'
    #segfile = 'segmentation_50_4.txt'
    #adj_file = 'Matthew.adjlist'
    #E_file = 'E_2.txt.npy'

    #data_dir = '../data/'
    #filename = 'Matthew_60min_sample.csv'
    #save_dir = '../result/OSC_Hurricane_Results/matthew/' + non + '/'
    #segfile = 'segmentation.txt'
    #adj_file = 'Matthew.adjlist'
    #E_file = 'E_5.txt.npy'

    #data_dir = '../data/'
    #filename = 'Irma_60min_sample.csv'
    #save_dir = '../result/OSC_Hurricane_Results/irma/' + non + '/'
    #segfile = 'segmentation.txt'
    #adj_file = 'Irma.adjlist'
    #E_file = 'E_1.txt.npy'

    data_dir = '../data/'
    filename = 'Harvey_60min_sample.csv'
    save_dir = '../result/Harvey_expl/'
    segfile = 'segV_lam1_2_lam2_0.1_lam3_0.1_clusV_3.csv'
    adj_file = 'Harvey.adjlist'
    E_file = 'E_0.2.txt'
    savefile = 'plot'
    clique = True
    thres = 0.8
    hurricane= 'Harvey'
    plot_result(data_dir, clique, adj_file, filename, E_file, segfile, save_dir, savefile, thres, hurricane)
