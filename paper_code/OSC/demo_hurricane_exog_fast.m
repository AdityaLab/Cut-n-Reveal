%  Hurricane demo fast version with exogenous inputs.

paths = genpath('libs/ncut');

addpath(paths);

%% Load Hurricane Data

rng(1);

%X = csvread('data/hurricane_matthew_normalized.csv');
X = csvread('data/hurricane_matthew_normalized_60min_sampling_interval.csv');
%X = csvread('data/hurricane_matthew_normalized_30min_sampling_interval.csv');
X=X(2:end,:); %Row 1 = Header row with column names. We ignore this.
X=X'; %Size(X) will be equal to 366 X 2478
%X = X/255; %no need for this since X is already column-normalized.

Y = csvread('data/hurricane_matthew_60min_sampling_interval.csv');
Y = Y(2:end,:);
Y=Y';
%X=Y;  %Normalized performs better than non-normalized.
corruption = 0;

w = randn(size(X)) * corruption;
X = X + w;

nbCluster = 4;

%% OSC

maxIterationbCluster = 300;
lambda_1 = .5;
lambda_2 = 100;  %check values around 1000.
gamma_1 = 0.01;
gamma_2 = 0.01;
p = 1.1;
diag = 1;

tic;
[Z, R, funVal] = OSCExog_fast(X,L, lambda_1, lambda_2, gamma_1, gamma_2, p, maxIterationbCluster, diag);
toc;

[oscclusters,osceigenvectors,osceigenvalues] = ncutW((abs(Z)+abs(Z')),nbCluster); %normalized cuts.

clusters = denseSeg(oscclusters,1);
plotClusters(clusters);


% plotting stuff

%clustering of non-normalized timeseries.
oldclusternum=-1;
segmentindices= [];
figure;
hold on
for i=1:size(Y,1)
    plot(Y(i,:)); %plot county timeseries.
end



%get segment boundaries
for i=1:size(Y,2)
    currclusternumber = clusters(i);
    if i==1
        oldclusternum = clusters(i);
    elseif(oldclusternum ~= currclusternumber)
        oldclusternum = currclusternumber;
        segmentindices = [i,segmentindices];
    end
end

segmentindices=sort(segmentindices); %sort in ascending order.
y1=get(gca,'ylim');
x1=get(gca,'xlim');
clustercolor=['b','y','r','k','g'];
alphaval=0.2;
for i=1:size(segmentindices,2)-1
    startidx=segmentindices(i);
    fprintf('Startidx %d and cluster %d\n',startidx,clusters(startidx));
    endidx=segmentindices(i+1);
    clusternum=clusters(startidx);
    color=clustercolor(clusternum);
    p=patch([startidx endidx endidx startidx],[1 1 y1(2) y1(2)],color);
    p.FaceAlpha=alphaval;
end

%populate the area before the first segment
endidx=segmentindices(1); %get the end idx of the first segment.
xlimit=x1(1); %get the beginning of the xlimit of the plot
clusternum=clusters(1); %get the cluster number of the first timestep.
color=clustercolor(clusternum)
p=patch([xlimit,endidx,endidx,xlimit],[1,1,y1(2),y1(2)],color);
p.FaceAlpha=alphaval;

%populate the end segment
endidx=segmentindices(end); %get the end idx of the last segment.
xlimit=x1(2);
clusternum=clusters(endidx);
color=clustercolor(clusternum)
p=patch([endidx,xlimit,xlimit,endidx],[1,1,y1(2),y1(2)],color);
p.FaceAlpha=alphaval;


hold off
%fileID = fopen('clusters_temp.csv','w');
%fprintf(fileID,'%d\n',clusters);
%fclose(fileID);

rmpath(paths);