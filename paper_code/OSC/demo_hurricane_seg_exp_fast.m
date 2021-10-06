%  Hurricane demo

paths = genpath('libs/ncut');
paths2 = genpath('libs/linspecer');

addpath(paths);
addpath(paths2);

%% Load Hurricane Data

rng(1);

inputdata='data/syn1_interp.csv'
plotting_input_data='data/syn1_interp.csv'
smatrixpath='HurricaneTimeseriesSegmentation/source/baseline/SubspaceClustering/data'

%Hurricane Data
%X = csvread('data/hurricane_matthew_normalized.csv',1);
%X = csvread('data/hurricane_matthew_normalized_60min_sampling_interval.csv',1);
X = csvread(inputdata,1);
%X = csvread('data/hurricane_matthew_60min_sampling_interval.csv');

%X=X(2:end,:); %Row 1 = Header row with column names. We ignore this.
X=X'; %Size(X) will be equal to 366 X 2478
%X = X/255; %no need for this since X is already column-normalized.

%Non Intrusive Load Monitoring Data. 
%X = csvread('data/non_intrusive_load_monitoring/non_intrusive_load_monitoring_dataset_1.csv',2);
%X = csvread('data/non_intrusive_load_monitoring/non_intrusive_load_monitoring_dataset_2.csv',2);
%X = csvread('data/non_intrusive_load_monitoring/non_intrusive_load_monitoring_dataset_3.csv',2);
%X = X';

%Normalization.
% means=mean(X,2);
% stdevs=std(X,0,2);
% for i=1:size(X,1)
%     X(i,:)= (X(i,:) - means(i))./stdevs(i);
%     
% end

% corruption = 0;
% 
% w = randn(size(X)) * corruption;
% X = X + w;

%for now we create a sample S matrix but in actuality it will be read in
%from a file.
 S = zeros(size(X,2),size(X,2)-1);
 S(:,100) = 1;
 %S(:,1:100:end) = 1; %every 50th column will be all ones.
%

% S matrix is to be provided as a csv file. 
%S is an m X m - 1  where m is the number of time steps. A column is
% assumed to be either all 1's or all 0's in the S matrix and for now we
% only consider integer elements.
if exist(strcat(smatrixpath,'s_matrix.txt'))==2
    S=csvread(strcat(smatrixpath,'s_matrix.txt')); 
end

nbCluster = 4;

%% OSC
maxIterationbCluster = 300;
lambda_1 = .5;
beta_1 = 0.000000005;  %making beta_1 low will erase the effect of the explanation.
beta_2 = 0.01;
p = 1.1;
diag = 1;

tic;
[Z, R, funVal] = OSC_seg_exp_fast(X,S,lambda_1,beta_1, beta_2, p, maxIterationbCluster, diag);
toc;

[oscclusters,osceigenvectors,osceigenvalues] = ncutW((abs(Z)+abs(Z')),nbCluster); %normalized cuts.

clusters = denseSeg(oscclusters,1);
%plotClusters(clusters);

%clustering of non-normalized timeseries.
oldclusternum=-1;
segmentindices= [];

%get segment boundaries
for i=1:size(X,2)
    currclusternumber = clusters(i);
    if i==1
        oldclusternum = clusters(i);
    elseif(oldclusternum ~= currclusternumber)
        oldclusternum = currclusternumber;
        segmentindices = [i,segmentindices];
    end
end

segmentindices=sort(segmentindices); %sort in ascending order.

%write segment indices to file for use by Explanation Algorithm. 
%Row Vector of segment indices ex: [4,50,100] indicates 4 segments 
% segment 1 : timestep 1 - 4 
% segment 2 : timestep 5 - 50
% segment 3 : timestep 51 - 100
% segment 4 : timestep 100 - end
csvwrite('temporal_segments.txt',segmentindices); 

%%% Start of  plotting Code. Can be commented out.
figure;
N=12;
tempvar = linspecer(N,'qualitative'); %get C for up to 12 different colors.
hold on
Y = csvread(plotting_input_data,1);
%Y = csvread('data/hurricane_matthew_60min_sampling_interval.csv');
Y = Y';
for i=1:size(Y,1)
    coloridx=mod(i,N);
    if coloridx==0
        coloridx=N;
    end
    plot(Y(i,:),'color',tempvar(coloridx,:),'Linewidth',1.6); %plot county timeseries.
end
xlim([0,size(Y,2)]);

y1=get(gca,'ylim');
x1=get(gca,'xlim');
clustercolor=['b','y','r','k','g'];
alphaval=0.00001;
for i=1:size(segmentindices,2)-1
    startidx=segmentindices(i);
    fprintf('Startidx %d and cluster %d\n',startidx,clusters(startidx));
    endidx=segmentindices(i+1);
    clusternum=mod(clusters(startidx),size(clustercolor,2));
    if clusternum==0
        clusternum=size(clustercolor,2);
    end
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
color=clustercolor(clusternum);
p=patch([endidx,xlimit,xlimit,endidx],[1,1,y1(2),y1(2)],color);
p.FaceAlpha=alphaval;

set([gca],'FontSize', 18);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold off
%fileID = fopen('clusters_temp.csv','w');
%fprintf(fileID,'%d\n',clusters);
%fclose(fileID);

rmpath(paths);
