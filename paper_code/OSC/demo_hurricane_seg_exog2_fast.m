%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_hurricane_seg_exog2_fast(resultsoutputdir,inputfile,visfile)
    %resultsoutputdir = '../result/Harvey_exog/';
    %inputfile = '../data/Harvey_60min_sample_normalized.csv';
    %visfile = '../data/Harvey_60min_sample.csv';
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
    xw= size(Y,2);
	%Hyperparameter.
    %lam2=[100,600,1200,1500];
    %lam2=[100000,120000,200000,500000];
    %nbclusters=[2,3,4];
    lam1=[2];
    lam2=[2];
    lam3=[2];
    maxiter=[300];
    nbclustersV=[4];
    %nbclustersU=[4];
    L = importdata('../Harvey_laplace.mat');
    [~,ln,~] = size(L);
	for iter=maxiter
		for nbCluster=nbclustersV
            for lambda_1=lam1 
			    for lambda_2=lam2
                            %% OSC
                            maxIterationbCluster = iter;
                            %lambda_1 = 0.5;
                            gamma_1 = 0.01;
                            gamma_2 = 0.01;
                            gamma_3 = 0.01; %for U matrix
                            beta_1 = 0.01; %for U matrix
                            p = 1.1;
                            diag = 1;
							tic;
                            [V, R, funVal1, iteration] = OSC_fast(X, lambda_1, lambda_2, gamma_1, gamma_2, p,... 
                                             maxIterationbCluster, diag);
                            toc;
						
                            [oscVclusters,oscVeigenvectors,oscVeigenvalues] = ncutW((abs(V)+abs(V')),nbCluster); %normalized V cuts.
                            clustersV = denseSeg(oscVclusters,1);
                            %disp('V cluster matrix');
                            %disp(clustersV);
                            %segments=sort(clustersV);
                            segmentindicesV= [];
                            oldclusternum=-1;
                            %get segment boundaries for V
                            for j=1:size(X,2)
                                currclusternumber = clustersV(j);%segments(j);
                                if j==1
                                    oldclusternum = clustersV(j);%segments(j);
                                elseif(oldclusternum ~= currclusternumber)
                                    oldclusternum = currclusternumber;
                                    segmentindicesV = [j,segmentindicesV];
                                end
                            end
                            segmentindicesV=sort(segmentindicesV); %sort in ascending order.
                            disp('segment indicesV');
                            disp(segmentindicesV);
                            %for each segment of V calculate U
                            startidx=1;
							endidx=xw;%1;
                            for lambda_3=lam3 
                                segmentindicesU= [];
								U = zeros(size(segmentindicesV,2)+1,ln);
                                for j=1:size(segmentindicesV,2)+1
                                    %startidx=endidx;
                                    %if j==size(segmentindicesV,2)+1
                                     %   endidx= xw;
                                    %else
                                    %    endidx= segmentindicesV(j);
                                    %end
                                    [U_idx, funValU, iterationU] = OSC_Exog2(X(:,startidx:endidx), L,...
                                        lambda_3, gamma_3, beta_1, p, iter, diag);
                                    [oscUclusters,oscUeigenvectors,oscUeigenvalues] = ncutW((abs(U_idx)+abs(U_idx')),nbCluster); 
                                                         %normalized U cuts
                                    clustersU = denseSeg(oscUclusters,1);
                                    fprintf('segment %d\n',j);
                                    U(j,:) = clustersU;
                                    %disp(U(j,:));
									segmentU= [];
									oldclusternum=-1;
									%get segment boundaries for U
									for k=1:size(X,1)
										currclusternumber = clustersU(k);%segments(j);
										if k==1
											oldclusternum = clustersU(k);%segments(j);
										elseif(oldclusternum ~= currclusternumber)
											oldclusternum = currclusternumber;
											segmentU = [k,segmentU];
										end
									end
									segmentU=sort(segmentU); %sort in ascending order.
									segmentindicesU= [j,segmentU];
									disp('segment indicesU');
									disp(size(segmentU));
									disp(segmentU);
                                end                               
                                %colorVar=['r','g','b','c'];
                                endidx=1;
                                figure;
                                hold on
                                xlim([0,size(Y,2)]);
                                for segidx=1:size(segmentindicesV,2)+1
                                    startidx=endidx;
                                    if segidx==(size(segmentindicesV,2))+1
                                        endidx= xw;
                                    else
                                        endidx= segmentindicesV(segidx);
                                    end
                                    %fprintf('segidx %d start %d end %d timesize %d\n',segidx, startidx, endidx,size(Y,2));
                                    N=200;
									tempvar = linspecer(N,'qualitative'); %get C for up to 12 different colors.
									color=0;
									oldclusternum=-1;
                                    for k=1:size(Y,1)
										currclusternumber = clustersU(k);
										if k==1
											oldclusternum = clustersU(k);
										elseif(oldclusternum ~= currclusternumber)
											oldclusternum = currclusternumber;
											color=color+1;
										end
										coloridx = mod(color,N);
										if coloridx==0
											coloridx=N;
                                        end
                                        plot(Y(k,startidx:endidx),'color',tempvar(coloridx,:)) %plot county timeseries.
                                    end
									disp('segment index U total');
									disp(color);
                                end
                                %xlim([0,size(Y,2)]);
                                resultfilename=strcat('lambda_1_',mat2str(lambda_1),'_lambda_2_',mat2str(lambda_2),...
                                '_lambda_3_',mat2str(lambda_3),'_numiter_',mat2str(maxIterationbCluster),...
                                '_numcluster_',mat2str(nbCluster)); 
                                
                                fignameV=strcat('osc_segments_V_U_',resultfilename);
                                figurefilenameV=strcat(resultsoutputdir,fignameV,'.pdf');
                
                                y1=get(gca,'ylim');
                                x1=get(gca,'xlim');
                                clustercolor=['b','y','r','m','g'];
                                %{
								alphaval = 0.01;
                                for j=1:size(segmentindicesV,2)-1
                                        startidx=segmentindicesV(j);
                                        endidx=segmentindicesV(j+1);
                                        clusternum=mod(clustersV(startidx),size(clustercolor,2));
                                        if clusternum==0
                                            clusternum=size(clustercolor,2);
                                        end
                                        color=clustercolor(clusternum);
                                        p=patch([startidx endidx endidx startidx],[1 1 y1(2) y1(2)],color);
                                        p.FaceAlpha=alphaval;
                                        line([segmentindicesV(j),segmentindicesV(j)],y1,'Color','k');
                                end
								%}
								%{
                                %populate the area before the first segment
                                 endidx=segmentindicesV(1); %get the end idx of the first segment.
                                 xlimit=x1(1); %get the beginning of the xlimit of the plot
                                 clusternum=clustersV(1); %get the cluster number of the first timestep.
                                 color=clustercolor(clusternum);
                                 p=patch([xlimit,endidx,endidx,xlimit],[1,1,y1(2),y1(2)],color);
                                 p.FaceAlpha=alphaval;

                                 %populate the end segment
                                 last = size(segmentindicesV,2);
                                 endidx=segmentindicesV(last); %get the end idx of the last segment.
                                 xlimit=x1(2);
                                 clusternum=clustersV(endidx);
                                 color=clustercolor(clusternum);
                                 p=patch([endidx,xlimit,xlimit,endidx],[1,1,y1(2),y1(2)],color);
                                 p.FaceAlpha=alphaval;
                                 line([segmentindicesV(last),segmentindicesV(last)],y1,'Color','k');
                                %}
								
                                for j=1:size(segmentindicesV,2)
                                    line([segmentindicesV(j),segmentindicesV(j)],y1,'Color','m');
                                end
                                
                                set([gca],'FontSize', 18);
                                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                                saveas(gcf,figurefilenameV);
                                print(figurefilenameV,'-dpdf','-bestfit');
                                hold off
                                                          
                                fileID = fopen('clusters_tempV.csv','w');
                                fprintf(fileID,'%d\n',clustersV);
                                fclose(fileID);
                                close(gcf);
                            
                                fileID = fopen('clusters_tempU.csv','w');
                                for j=1:size(segmentindicesV)+1
                                    fprintf(fileID,'%d\n',U(j,:));
                                    fprintf(fileID,'\n');
                                end
                                fclose(fileID);
                            end
                            
                             

                            
                            %segment
                            segoutputfileV=strcat('osc_segment_indices_V_',resultfilename);
                            %segoutputfileU=strcat('osc_segment_indices_U_',resultfilename);
                            segmentfilenameV=strcat(resultsoutputdir,segoutputfileV,'.csv');
                            %segmentfilenameU=strcat(resultsoutputdir,segoutputfileU,'.csv');
                            csvwrite(segmentfilenameV,segmentindicesV);
                            %csvwrite(segmentfilenameU,segmentindicesU);
                            
                            %V matrix write
                            
                            vmatrixoutputfile=strcat(resultsoutputdir,'V_matrix_',resultfilename,'.csv');
                            fprintf('vmatrixoutputfile %s\n',vmatrixoutputfile);
                            dlmwrite(vmatrixoutputfile,V);
                            %{
                            %U matrix write
                            umatrixoutputfile=strcat(resultsoutputdir,'U_matrix_',resultfilename,'.csv');
                            fprintf('umatrixoutputfile %s\n',umatrixoutputfile);
                            dlmwrite(umatrixoutputfile,U);

                            %Affinity matrix V write
                            
                            affinitymatrixoutputfileV=strcat(resultsoutputdir,'V_affinity_matrix_',resultfilename,'.csv');
                            dlmwrite(affinitymatrixoutputfileV,(abs(V)+abs(V')));
                           
                            %Affinity matrix U write
                            affinitymatrixoutputfileU=strcat(resultsoutputdir,'U_affinity_matrix_',resultfilename,'.csv');
                            dlmwrite(affinitymatrixoutputfileU,(abs(U)+abs(U')));

                            %figure output name V
                            %}
                end
		    end
        end
    end
	rmpath(paths);
	rmpath(paths2);
end
