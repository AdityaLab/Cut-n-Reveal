#hyperparameters to obtain segmentation by running matlab file
lam1=[0.1]
lam2=[0.1]
lam3 =[0.1]
beta=[0.7 1]
nbclustersU=[2]
nbclustersV=[3]
K='3'
segLimit='10'
MATLAB_path = '/usr/local/bin/./matlab'
Result_file='../result/'
Data_file='../data/Harvey_60min_sample_normalized.csv'
Vis_file='../data/Harvey_60min_sample.csv'
Laplace_file='../data/Harvey_laplace.mat'

demo:
	$(MATLAB_path) -nodisplay -nosplash -nodesktop -r "demo_spatial_seg_exog3_fast($(Result_file),$(Data_file),$(Vis_file),$(Laplace_file),$(K),$(segLimit),$(lam1),$(lam2),$(lam3),$(beta),$(nbclustersV),$(nbclustersU))"

	