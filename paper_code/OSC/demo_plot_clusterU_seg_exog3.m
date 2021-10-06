% find clusters of U from UU'
function [] = demo_plot_clusterU_seg_exog3(resultsoutputdir,inputfile,resultfilename)
    paths = genpath('libs/ncut');
	paths2 = genpath('libs/linspecer');
    fprintf('results dir %s\n',resultsoutputdir);
	addpath(paths);
	addpath(paths2);
nbClustersU =[2,3,4];
for nbClusterU=nbClustersU
    affineU = csvread(inputfile);
    disp(size(affineU));
    [oscUclusters,oscUeigenvectors,oscUeigenvalues] = ncutW(affineU,nbClusterU); %normalized U cuts
    clustersU = denseSeg(oscUclusters,1);
    clusterfilenameU=strcat(resultsoutputdir,'clusters_U_',mat2str(nbClusterU),'_',resultfilename,'.csv');
    fileID = fopen(clusterfilenameU,'w');
    fprintf(fileID,'%d\n',clustersU);
    fclose(fileID);
end

rmpath(paths);
rmpath(paths2);
end