%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_hurricane_fast_func(resultsoutputdir,inputfile,visfile)
	paths = genpath('libs/ncut');
	paths2 = genpath('libs/linspecer');
    fprintf('results dir %s',resultsoutputdir);
	addpath(paths);
	addpath(paths2);

	%% Load Hurricane Data

	rng(1);

	X = csvread(inputfile,1);
	Y = csvread(visfile,1); 

	X=X'; %Size(X) will be equal to 366 X 625
	Y=Y'; %This is the data to visualize 366 X 625.
	corruption = 0;

	w = randn(size(X)) * corruption;
	X = X + w;

	%Hyperparameter.
        %lam2=[100,600,1200,1500];
        %lam2=[100000,120000,200000,500000];
        %lam1 =[0.5];
        %lam2=[0.5];
        lam1=[0.1,0.5,2];
        lam2=[0.1,0.5,2];
        maxiter=[300];
        nbclusters=[2,3,4];
        kernels={'original'};
    fileID = fopen('loss.txt','w');
    fprintf(fileID,'lambda1, lambda2, gam1, gam2, cluster, loss');
	for iter=maxiter
		for nbCluster=nbclusters
            for lambda_1=lam1 
			    for lambda_2=lam2
                    for i=1:size(kernels,2)
                    %% OSC
                    maxIterationbCluster=iter
                    %lambda_1 = 0.5;
                    gamma_1 = [0.03,.05,.07];
                    gamma_2 = [0.03,.05,.07];
                    p = 1.1;
                    diag = 1;
                    kernel=kernels(i);
                    for gam1=gamma_1
                        for gam2=gamma_2
                            tic;
                            [Z, R, funVal,it] = OSC_fast(X, lambda_1, lambda_2, gam1, gam2, p, maxIterationbCluster, diag);
                            toc;
                            
                            fprintf(fileID,'%.2f,%.2f,%.2f,%.2f,%d,%.3f\n',lambda_1, lambda_2, gam1, gam2, nbCluster, funVal(it));
                    
                            if strcmp(kernel,'original')
                                [oscclusters,osceigenvectors,osceigenvalues] = ncutW((abs(Z)+abs(Z')),nbCluster); %normalized cuts.
                            elseif strcmp(kernel,'multiplication')
                                [oscclusters,osceigenvectors,osceigenvalues] = ncutW(Z*Z',nbCluster); %normalized cuts.
                            end
                            clusters = denseSeg(oscclusters,1);
                            %plotClusters(clusters);
                            % nikhil plotting stuff
                            %clustering of non-normalized timeseries.
                            oldclusternum=-1;
                            segmentindices= [];
                            
                            %get segment boundaries
                            for j=1:size(X,2)
                                currclusternumber = clusters(j);
                                if j==1
                                    oldclusternum = clusters(j);
                                elseif(oldclusternum ~= currclusternumber)
                                    oldclusternum = currclusternumber;
                                    segmentindices = [j,segmentindices];
                                end
                            end

                            segmentindices=sort(segmentindices); %sort in ascending order.
                            fV = prune_candidate_segment(.05,segmentindicesV,nbClusterV,size(X,2),1);
                            if fV
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
                                %disp(segmentindices)
                                resultfilename=strcat('lambda_1_',mat2str(lambda_1),'_lambda_2_',mat2str(lambda_2),'_gamma_1_',mat2str(gam1),'_gamma_2_',mat2str(gam2),'_numcluster_',mat2str(nbCluster));  
                                segoutputfile=strcat('osc_segment_indices_',resultfilename);     
                                segmentfilename=strcat(resultsoutputdir,segoutputfile,'.csv');
                                csvwrite(segmentfilename,segmentindices);
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
                                    p1=patch([startidx endidx endidx startidx],[1 1 y1(2) y1(2)],color);
                                    p1.FaceAlpha=alphaval;
                                end
                                %populate the area before the first segment
                                endidx=segmentindices(1); %get the end idx of the first segment.
                                xlimit=x1(1); %get the beginning of the xlimit of the plot
                                clusternum=clusters(1); %get the cluster number of the first timestep.
                                color=clustercolor(clusternum)
                                p1=patch([xlimit,endidx,endidx,xlimit],[1,1,y1(2),y1(2)],color);
                                p1.FaceAlpha=alphaval;
                                %populate the end segment
                                endidx=segmentindices(end); %get the end idx of the last segment.
                                xlimit=x1(2);
                                clusternum=clusters(endidx);
                                color=clustercolor(clusternum);
                                p1=patch([endidx,xlimit,xlimit,endidx],[1,1,y1(2),y1(2)],color);
                                p1.FaceAlpha=alphaval;
                                set([gca],'FontSize', 18);
                                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                                %saveas(gcf,figurefilename);
                                print(figurefilename,'-dpdf','-bestfit')
                                hold off
                                %fileID = fopen('clusters_temp.csv','w');
                                %fprintf(fileID,'%d\n',clusters);
                                %fclose(fileID);
                                close(gcf);
                            end
                        end
                    end
                  end
			    end
		      end
             end
           end
	rmpath(paths);
	rmpath(paths2);
    fclose(fileID);
end
