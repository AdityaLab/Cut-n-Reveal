%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_hurricane_clustering_and_visualization(resultsoutputdir,lambda1,lambda2,numcluster,visfile,kernel)
    paths = genpath('libs/ncut');
    paths2 = genpath('libs/linspecer');
    fprintf('results dir %s',resultsoutputdir);
    addpath(paths);
    addpath(paths2);

    %% Load Hurricane Data

    rng(1);
    
    
    Y = csvread(visfile,1); 

   
    Y=Y'; %This is the data to visualize 366 X 625.


    %% OSC
    maxIterationbCluster=300;
    nbCluster=numcluster;
    lambda_1 = lambda1;
    lambda_2 = lambda2;
    gamma_1 = 0.01;
    gamma_2 = 0.01;
    p = 1.1;
    diag = 1;
    
    vmatrixinputfilename=strcat('lambda_1_',mat2str(lambda_1),'_lambda_2_',mat2str(lambda_2),'_numiter_',mat2str(maxIterationbCluster),'_numcluster_',mat2str(nbCluster));  
    vmatrixinputfile=strcat(resultsoutputdir,'V_matrix_',vmatrixinputfilename,'.csv');
    Z = csvread(vmatrixinputfile);
    if kernel=='original'
        [oscclusters,osceigenvectors,osceigenvalues] = ncutW((abs(Z)+abs(Z')),nbCluster); %normalized cuts.
    elseif kernel=='maxtirx_multiplication'
        [oscclusters,osceigenvectors,osceigenvalues] = ncutW(Z*Z',nbCluster); %normalized cuts.
    end 
    clusters = denseSeg(oscclusters,1);


    % plotting stuff
    %clustering of non-normalized timeseries.
    oldclusternum=-1;
    segmentindices= [];
     figure;
     N=12;
     tempvar = linspecer(N,'qualitative'); %get C for up to 12 different colors.
     hold on

     for i=1:size(Y,1)
         coloridx=mod(i,N);
         if coloridx==0
         coloridx=N;
         end
         plot(Y(i,:),'color',tempvar(coloridx,:),'Linewidth',1.6); %plot county timeseries.
     end
     xlim([0,size(Y,2)]);


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
    resultfilename=strcat('_affinitymatrix_',kernel,'_lambda_1_',mat2str(lambda_1),'_lambda_2_',mat2str(lambda_2),'_numiter_',mat2str(maxIterationbCluster),'_numcluster_',mat2str(nbCluster));  

    %segment
    segoutputfile=strcat('osc_segment_indices_',resultfilename);     
    segmentfilename=strcat(resultsoutputdir,segoutputfile,'.csv');
    csvwrite(segmentfilename,segmentindices);

    %V matrix write
%     vmatrixoutputfile=strcat(resultsoutputdir,'V_matrix_',resultfilename,'.csv');
%     fprintf('vmatrixoutputfile %s\n',vmatrixoutputfile);
%     dlmwrite(vmatrixoutputfile,Z);

    %Affinity matrix write
    affinitymatrixoutputfile=strcat(resultsoutputdir,'V_affinity_matrix_',resultfilename,'.csv');
    
    dlmwrite(affinitymatrixoutputfile,(abs(Z)+abs(Z')));

     %figure output name
     figname=strcat('osc_segments_',resultfilename);
     figurefilename=strcat(resultsoutputdir,figname,'.pdf');


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
     %saveas(gcf,figurefilename);
     print(figurefilename,'-dpdf','-bestfit')
     hold off
     
     close(gcf);

    rmpath(paths);
    rmpath(paths2);
end
