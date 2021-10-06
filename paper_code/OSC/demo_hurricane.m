%  Faces demo

paths = genpath('libs/ncut');

addpath(paths);

%% Load faces

rng(1);

%X = csvread('data/hurricane_matthew_normalized.csv');
X = csvread('data/hurricane_matthew_normalized_60min_sampling_interval.csv');
X=X'; %Size(X) will be equal to 366 X 2478
%X = X/255; %no need for this since X is already column-normalized.

corruption = 0;

w = randn(size(X)) * corruption;
X = X + w;

nbCluster = 2;

%% OSC

maxIterationbCluster = 300;
lambda_1 = .5;
lambda_2 = 5;
gamma_1 = 0.01;
gamma_2 = 0.01;
p = 1.1;
diag = 1;

tic;
[Z, R, funVal] = OSC(X, lambda_1, lambda_2, gamma_1, gamma_2, p, maxIterationbCluster, diag);
toc;

[oscclusters,~,~] = ncutW((abs(Z)+abs(Z')),nbCluster);



clusters = denseSeg(oscclusters,1);
plotClusters(clusters);

rmpath(paths);