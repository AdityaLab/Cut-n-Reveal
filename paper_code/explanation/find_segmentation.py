import os
import sys
import random
import numpy as np
import networkx as nx
import subprocess
import scipy.io
import time
matlab = '/raid/R2014A/bin/matlab'

def find_segmentation(data_dir, hurricane, k, alpha, lmda, 
    save_dir, step, thres, nclusters, max_iter_seg, max_iter):
    """
    the algorithm to find segmentation and the explanation in an iterative manner, the explanation part is solved using a closed form solution (not using a closed form solution becuase it has negative value, now using a quadratic programming solver)
    One can easily change the code to select how to find the segmentation in each iteration, either using a naive algorithm or using Nikhil's formular
    k: number of cut points 
    alpha: the weight for the laplacian regularization
    lmda: the weight for the lasso regularization
    step: the minimum number of time units between two adjacent cut points (like the resolution of the segmentation)
    thres: threshold for the segmentation changes, if the change is less than the threshold, stop the iteration process
    nclusters: the number of clusters to use in the spectral clustering algorithm to find the segmentation
    max_iter_seg: the maximum number of itereation for finding the segmentation
    max_iter: the maximum number of total iterations between finding the segmentation and the explanation
    """
    k = int(k)
    save_dir += hurricane + '/' + str(k) + '/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    print 'reading the data'
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
    for i in range(len(lines)):
        line = lines[i]
        data.append([float(x) for x in line.strip().split(',')])
    data = np.array(data)
    l = len(lines)
    print 'reading the network'
    G = nx.read_adjlist(data_dir + hurricane + '.adjlist', delimiter = '\t')
    L = nx.laplacian_matrix(G, counties).toarray()

    print 'init segmentation, and difference'
    pre_segmentation = init_seg(l, k, step)
    diff = get_diff(data, l, k, pre_segmentation)
    sf = open(save_dir + 'seg_changes_iterations', 'wb')
    for i in range(max_iter):
        print 'iteration ' + str(i)
        E = callMatlab(diff, L, alpha, lmda, save_dir)
        #S = new_seg(data, pre_segmentation, E, step)
        S = new_seg_2(data, pre_segmentation, E, step, save_dir, data_dir, hurricane)
        diff = get_diff(data, l, k, S)
        S_diff = seg_diff(data, pre_segmentation, S) 
        print S_diff
        sf.write(str(i) + '\t' + str(S_diff) + '\n')
        if S_diff< thres:
            break
        pre_segmentation = S
    sf.close()
    print 'saving results'
    save_result(S, E, k, save_dir)


def callMatlab(diff, L, alpha, lmda, save_dir):
    #D = '['
    #for d in diff:
    #    D += ','.join([str(x) for x in d]) + ';'
    #D = D[:-1] + '];'
    #LL = '['
    #for l in L:
    #    LL += ','.join([str(x) for x in l]) + ';'
    #LL = LL[:-1] + '];'
    #matlab_cmd = 'diff = ' + D
    #matlab_cmd += 'L = ' + LL
    #matlab_cmd += 'find_exp(diff, L, ' + str(alpha) + ',' + str(lmda) + ',\'' + save_dir + '\');exit;'
    scipy.io.savemat(save_dir + 'dif.mat', mdict = {'dif': diff})
    scipy.io.savemat(save_dir + 'L.mat', mdict = {'L': L})
    matlab_cmd = 'find_exp(' + str(alpha) + ',' + str(lmda) + ',\'' + save_dir + '\');exit;'
    print("matlab -nosplash -nojvm -r ",matlab_cmd)
    subprocess.call([matlab, "-nosplash", "-nojvm", "-r", matlab_cmd])
    with open(save_dir + 'E.txt') as f:
        lines = f.readlines()
    E = []
    for line in lines:
        E.append([float(x) for x in line.strip().split('\t')])
    E = np.array(E)
    #normalize E
    #E = E / E.sum(axis = 0)
    #E = (E - E.min(0)) / (E.max(0) - E.min(0))
    return E
def new_seg_2(data, pre_S, E, step, save_dir, data_dir, hurricane):
    Sexp = []
    for i in range(len(data) - 1):
        on_cut, cur_cut, prev_cut, next_cut = locate(i + 1, pre_S)
        if on_cut:
            exp = E[:, cur_cut]
            Sexp.append(time_stamp_diff(data[i + 1], data[i], exp))
        elif prev_cut == -1:
            exp = E[:, 0]
            Sexp.append(time_stamp_diff(data[i + 1], data[i], exp))
        elif next_cut == len(pre_S):
            exp = E[:, len(pre_S) - 1]
            Sexp.append(time_stamp_diff(data[i + 1], data[i], exp))
        else:
            prev_exp = E[:, prev_cut]
            next_exp = E[:, next_cut]
            prev_dif = time_stamp_diff(data[i + 1], data[i], prev_exp)
            next_dif = time_stamp_diff(data[i + 1], data[i], next_exp)
            if abs(prev_dif) > abs(next_dif):
                Sexp.append(prev_dif)
            else:
                Sexp.append(next_dif)
    #Sexp = Sexp > np.mean(Sexp)
    Sexp = mask(Sexp)
    Sexp = [1 if x else 0 for x in Sexp]
    #sf = open(save_dir + 's_matrix.txt', 'wb')
    sf = open(save_dir + 's_matrix.txt', 'wb')
    for i in range(len(data)):
        sf.write(','.join([str(x) for x in Sexp]) + '\n')
    sf.close()
    #matlab_cmd = 'demo_hurricane_seg_exp_fast;exit;'
    matlab_cmd = 'seg_exp_fast(\'' + data_dir + hurricane + '_interp.csv\','
    matlab_cmd += '\'' + save_dir + '\', ' + str(nclusters) + ',' + str(max_iter_seg) + ');exit;' 
    print matlab_cmd
    subprocess.call([matlab, "-nosplash", "-nojvm", "-r", matlab_cmd])
    with open(save_dir + 'temporal_segments.txt')as f:
        line = f.readlines()[0]
    new_S = line.strip().split(',')
    new_S = [int(x) for x in new_S]
    return new_S

def mask(s):
    """
    keep the top 15% values and make them 1, the others 0
    """
    x = zip(range(len(s)), s)
    x = sorted(x, key = lambda a:a[1], reverse = True)
    thres = sum(s) * 0.15
    su = 0.0
    result = [0] * len(s)
    for ind, val in x:
        s[ind] = 1
        su += val
        if su > thres:
            break
    return result
def time_stamp_diff(d1, d2, exp):
    """
    Calculate the difference of the data in two adjacent time stamp based on the explanation.
    d1 is the later time stamp, d2 is for the previous time stamp
    """
    s = 0.0
    for i in range(len(d1)):
        s += (d1[i] - d2[i]) * exp[i]
    return s
def locate(t, S):
    """
    find the location of time stamp i in a sgementation S.
    if t is on a cut point, return on_cut True
    otherwise return the prev cut point and the next cut point that surround t
    """
    for i in range(len(S)):
        s = S[i]
        if t == s:
            return True, i, None, None
        if t < s:
            return False, None, i - 1, i
    return False, None, i, i + 1
def new_seg(data, pre_S, E, step):
    S = []
    for i in range(len(pre_S)):
        if i == 0:
            start = 1 
        else:
            start = S[-1] + 1
        if i == len(pre_S) - 1:
            end = len(data) - 1
        else:
            end = pre_S[i + 1]
        #print start, end
        cand = range(start, end, step)
        #print cand
        max_dist = -sys.maxint
        max_c = None
        for c in cand:
            if start == 1:
                d1 = data[:c + 1]
            else:
                d1 = data[start: c + 1]
            d2 = data[c + 1: end + 1]
            dist = get_distance(d1, d2)
            E_dist = np.inner(E[:, i], dist)
            if E_dist > max_dist:
                max_dist = E_dist
                max_c = c
        S.append(max_c)
    return S

def seg_diff(data, S1, S2):
    """
    calculate the difference of two segmentations.
    """
    """
    if len(S1) != len(S2):
        raise 'the length of the two segments are not the same'
    d = 0.0
    for i in range(len(S1)):
        d += abs(S1[i] - S2[i])
    return d
    """
    s1 = 0.0
    for i in range(len(S1)):
        x = S1[i]
        min_d = sys.maxint 
        for j in range(len(S2)):
            y = S2[j]
            d = abs(x - y)
            if d < min_d:
                min_d = d
        s1 += min_d
    s1 /= len(S1)
    s2 = 0.0
    for i in range(len(S2)):
        x = S2[i]
        min_d = sys.maxint 
        for j in range(len(S1)):
            y = S1[j]
            d = abs(x - y)
            if d < min_d:
                min_d = d
        s2 += min_d
    s2 /= len(S2)
    return (s1 + s2) / 2
def save_result(S, E, k, save_dir):
    sf = open(save_dir + 'S_' + str(k) + '.txt', 'wb')
    sf.write(str(S))
    sf.close()
    np.save(save_dir + 'E_' + str(k) + '.txt', E)
    #sf = open(save_dir + 'E_' + str(k) + '.txt', 'wb')
    #sf.write(np.array_str(E))
    #sf.close()

def get_distance(d1, d2):
    """
    Given the set of time series in two time segments, calculate the distance in each time series.
    """
    noC = len(d1[0])#number of columns
    #dist_matrix = np.array([])
    dm_mean = []
    dm_std = []
    dm_grad = []
    for i in range(noC):
        ts1 = d1[:, i]
        ts2 = d2[:, i]
        #temp = [] 
        m1 = np.mean(ts1)
        m2 = np.mean(ts2)
        dm_mean.append(abs(m1 - m2))
        s1 = np.std(ts1)
        s2 = np.std(ts2)
        dm_std.append(abs(s1 - s2))
        if len(ts1) > 10:
            g1 = np.mean(np.gradient(ts1))
        else:
            g1 = 0
        if len(ts2) > 10:
            g2 = np.mean(np.gradient(ts2))
        else: 
            g2 = 0
        dm_grad.append(abs(g1 - g2))
        #dist_matrix = np.append(dist_matrix,temp)
    #if not max(dm_mean) - min(dm_mean) < 0.001:
    #    dm_mean = [(x - min(dm_mean)) * 1.0 /(max(dm_mean) - min(dm_mean)) for x in dm_mean]
    #dm_mean = [x * 1.0 /max(dm_mean) for x in dm_mean]
    #dm_std = [(x - min(dm_std)) * 1.0 /(max(dm_std) - min(dm_std)) for x in dm_std]
    #dm_std = [x * 1.0 /max(dm_std) for x in dm_std]
    #if not max(dm_grad) - min(dm_grad) < 0.001:
    #    dm_grad = [(x - min(dm_grad)) * 1.0 /(max(dm_grad) - min(dm_grad)) for x in dm_grad]
    #dm_grad = [x * 1.0 /max(dm_grad) for x in dm_grad]
    if len(d1) < 20 or len(d2) < 20:
        l_pen = 0.0
    else:
        l_pen = 1
    #return [l_pen * 0.5 * (dm_mean[i] + dm_std[i]) for i in range(len(dm_mean))]
    #return [l_pen * 0.5 * (dm_mean[i] + dm_grad[i]) for i in range(len(dm_mean))]
    #return [l_pen * dm_mean[i] for i in range(len(dm_mean))]
    return [l_pen * dm_grad[i] for i in range(len(dm_grad))]
    #dist_matrix = (dist_matrix - dist_matrix.min(0)) * 1.0 / dist_matrix.ptp(0)
    #return [np.mean(x) for x in dist_matrix]

def get_diff(data, l, k, s):
    """
    for each cut point in s, find the distance of the adjacent time segments.
    """
    diff = []
    k = len(s)
    for i in range(k):
        if i == 0:
            d1 = data[0:s[0] + 1]
            if k > 1:
                d2 = data[s[0] + 1: s[1] + 1]
            else:
                d2 = data[s[0] + 1:]
        elif i == k - 1:
            d1 = data[s[i - 1] + 1: s[i] + 1]
            d2 = data[s[i] + 1:]
        else:
            d1 = data[s[i - 1] + 1: s[i] + 1]
            d2 = data[s[i] + 1: s[i + 1] + 1]
        diff.append(get_distance(d1, d2))
    return np.array(diff)
        
def init_seg(l, k, step):
    """
    the returned set of cut points [c1, c2, c3, ...] represents the time segments of [0, c1], (c1, c2], (c2, c3], ... (ck, l - 1]  
    """
    if k > (l - 1) * 1.0 / step:
        raise 'k is too large'
    if k < 1:
        raise 'k is too small'
    poten = range(1, l - 1, step)
    random.shuffle(poten)
    return sorted(poten[:k])

if __name__ == '__main__':
    data_dir = './data/'
    #hurricane = 'Matthew'
    #hurricane = 'syn2'
    hurricane = sys.argv[1] #data input file
    k = sys.argv[2]  #num cutpoints we want, not very useful in combined version.
    alpha = 0.2 #liangzhe's part
    lmda = 0.001 #liangzhe's part
    save_dir = './result/'
    step = 20  #liangzhe's part
    thres = 40  #liangzhe's part
    nclusters = 4  #nikhil part
    max_iter_seg = 200 #nikhil part
    max_iter = 200 #combined formulation.
    find_segmentation(data_dir, hurricane, k, alpha, lmda, save_dir, step, thres, nclusters, max_iter_seg, max_iter)
