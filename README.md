'Cut-n-Reveal: Time-series Segmentations with Explanations (CnR-UV)', Nikhil Muralidhar, Liangzhe Chen, Anika Tabassum, Supriya Chinthavali, Naren Ramakrishnan, and B. Aditya Prakash
==========================================================================
Paper:
-------------
Both paper and supplementary readings are added here. Check paper.pdf and appendix.pdf

Usage:
Note: You need to set the correct MATLAB_path in the makefile (Including the MATLAB executable).
```
- Example:
    MATLAB_path = '/usr/local/bin/./matlab'
```
To run cnrUV for segmentation do as follows,
```
>> make demo  
```
'make demo' will run for the sample data (Harvey) given in data/ directory.

segmentation/: 
-------------
For getting segmentation. You can directly run the code using following command
```
- Example: 
    MATLAB_path -nodisplay -nosplash -nodesktop -r "demo_spatial_seg_exog3_fast(<result_file>,<data_file>,<vis_file>,<laplace_file>,k,segLimit,lam1,lam2,lam3,beta,nbclustersV,nbclustersU)"
```

Input: 
------
-- <data_file> .csv file located in data/ is a comma seperated file. Each row represents a multivariate observation. Each column represents a time-series. Total N columns (#time-series) and T rows (#timestamps). For input to get segmentation this data file needs to be normalized. First row should contain name of time-series.

The data.csv file should look like this:
```
ts1,ts2,ts3,...,tsn
x11,x21,x31,...,xn1
x12,x22,x32,...,xn2
.
.
.
x1T,x2T,x3T,...,xnT
```

--<laplace_file> .mat file as a laplacian matrix should be located in data/. If not present then generate this file using following command:
python laplacian.py <data_file> <adj_list_file> <save_file>  

	- <data_file> : Directory of the dataset
	- <adj_list_file>: Directory of the adjacency list file. Each row first column contains the time-series/region and next all columns contain its neghboring regions.
	- <save_file>: Directory for the results to be saved   

-- <vis_file> .csv located in data/ is a comma separated file. This can be same as data file or a non-normalized version of the data file. This file is to plot the segmentation.

-- <result_file> directory where results should be saved.

-- hyperparameters:
   
   -k: latent parameter. An integer number in string format. e.g., '5'.  
   
   -seg_limit: limit the max number of segments that can be considered for result. An integer number in string format, e.g., '10'.
   
   -lam1: an array of float values. first try between (0,1), e.g., [0.1] 
   
   -lam2: an array of float values. first try between (0,1), e.g., [0.1] 
   
   -lam3: an array of float values. first try between (0,1), e.g., [0.1] 
   
   -beta: an array of float values. first try between (0.5,1), e.g., [0.7] 
   
   -nbclustersV: an array of integer values. Represents the max number of temporal clusters, e.g., [3] 
   
   -nbclustersU: an array of integer values. Represents the max number of spatial clusters, e.g., [2]

Note: 
-A set of one hyperparameter should give desired result. To make things simpler, we select the candidate whose segmentation is closest to ground-truth or loss value in ../result/loss.txt is minimum (when no ground-truth).
- There are also other hyperparameters in the matlab file, gamma1, gamma2, gamma3 but they are robust and we set them to be 0.5 in all cases. They can also be changed if necessary.

Output:
-------
segV_lam1_v1_lam2_v2_lam3_v3_clusV_v4_l_v5_clusU_v6.csv: in the <result_file> directory. The main result file which gives the segmentation. The values in the file v1,v2,v3,v4,v5,v6 are the variable values of the associated set of hyperparameters. l is the value of k.

loss.txt: log the loss obtained for each set of hyperparameter values. This is necessary to select the desired segmentation. (We select the ones with minimum loss).

oscV_lam1_v1_lam2_v2_lam3_v3_clusV_v4_clusU_v5.pdf: The time-series plot with the segmentation obtained by associated hyperparameter.

U_affinity_matrix_lam1_v1_lam2_v2_lam3_v3_clusV_v4_l_3_clusU_v5.csv: U'U, where U is learnt by cnrUV. This is needed to obtain explanation.

clustersU_lam1_v1_lam2_v2_lam3_v3_clusV_v4_clusU_v5.csv: each row i represents which cluster time-series i should belong to. This just shows how time-series are clustered spatially.

explanation/: 
------------
To run cnrUV for explanation (which time-series are culprit for the segmentation) do as follows,
```
>> make demo 
```
Note:
--'make demo' will run for the same sample data (Harvey) given in data/ directory
--'make demo' should must run after running segmentation code. If no segmentation is found, this will not work

```
- Example: 
    python find_exp_new.py <matlab_path> <data_dir> <data_file> <result_file> <adj_list_file> <seg_file> <file_affinityU> alpha lmda nclusters
```

Input:
-------
--<data_dir>: data directory

--<data_file>: data file name

--<result_file>: directory to save

--<seg_file> : seg file name in the result directory (obtained from segmentation). This segmentation will be used for explanation.

--<file_affinity_U>: U_affinity_matrix file name found by running segmentation code. This output file will be same found from the corresponding <seg_file> (same hyperparams)

-- hyperparams:
    
    - alpha: a float value select how many to choose for explanation
    
    - lmbda: a float value
    
-- nclusters: number of cluster chosen for V (during segmentation)

Output:
-------
 #located in result_file directory

-E_alpha.txt: for each cutpoint explanation weight of every time-series

-plot_pdf: plot of time-series and segmentation

-plot_cut_impC.txt: txt file showing which time-series chosen for explanation at each cutpoint

-plot_exp: plot of time-series and segmentation+explanation

-plot_exp_i: plot of time-series and explanation of cutpoint i
