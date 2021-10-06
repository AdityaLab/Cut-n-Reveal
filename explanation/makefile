matlab='/usr/local/bin/./matlab'
data_dir = '../data/'
filename = 'Harvey_60min_sample.csv'
save_dir = '../result/'
adj_file = 'Harvey.adjlist'
segfile = 'segV_lam1_0.1_lam2_0.1_lam3_0.1_clusV_3_l_5.csv'
file_affinityU='U_affinity_matrix_lam1_0.1_lam2_0.1_lam3_0.1_clusV_3_l_3_clusU_3.csv'

alpha = 0.2 #liangzhe's part
lmda = 0.001 #liangzhe's part
nclusters = 3  #nikhil part, same values as nclusterV

demo:
	python find_exp_new.py $(matlab) $(data_dir) $(filename) $(save_dir) $(adj_file) $(segfile) $(file_affinityU) $(alpha) $(lmda) $(nclusters)