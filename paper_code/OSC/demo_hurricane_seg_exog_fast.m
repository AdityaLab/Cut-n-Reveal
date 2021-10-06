%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_hurricane_seg_exog_fast(resultsoutputdir,inputfile,visfile,laplaceFile)
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

	%Hyperparameter.
        lam1=[100,600,1200,1500];
        lam2=[100000,120000,200000,500000];
        %nbclusters=[2,3,4];
        %lam1=[100];
        %lam2=[200000];
        lam3=[.5,50,100];
        maxiter=[300];
        nbclusters=[4];
	for iter=maxiter
		for nbCluster=nbclusters
            for lambda_1=lam1 
			    for lambda_2=lam2
                    for lambda_3=lam3
                            %% OSC
                            maxIterationbCluster = iter;
                            %lambda_1 = 0.5;
                            gamma_1 = 0.5;
                            gamma_2 = 0.5;
                            gamma_3 = 0.5;
                            beta_1 = 0.01;
                            p = 1.1;
                            diag = 1;
                            L = importdata(laplaceFile);
                            tic;
                            [V, U, funVal] = OSCU_fast(X, L, lambda_1, lambda_2, lambda_3, gamma_1, gamma_2, gamma_3, p,... 
                                             beta_1, maxIterationbCluster, diag, diag);
                            toc;

                            [oscVclusters,oscVeigenvectors,oscVeigenvalues] = ncutW((abs(V)+abs(V')),nbCluster); %normalized V cuts.
                            [oscUclusters,oscUeigenvectors,oscUeigenvalues] = ncutW((abs(U)+abs(U')),nbCluster); %normalized U cuts
                            clustersV = denseSeg(oscVclusters,1);
                            %disp('V cluster matrix');
                            %disp(clustersV);
                            clustersU = denseSeg(oscUclusters,1);
                            %disp('U cluster matrix');
                            %disp(clustersU);
                            resultfilename=strcat('lambda_1_',mat2str(lambda_1),'_lambda_2_',mat2str(lambda_2),'_numiter_',mat2str(maxIterationbCluster),'_numcluster_',mat2str(nbCluster));  

                            % nikhil plotting stuff
                
                            %clustering of non-normalized timeseries.
                            oldclusternum=-1;
                            figure;
                            N=250;
                            tempvar = linspecer(N,'qualitative'); %get C for up to 12 different colors.
                            %colorVar=['r','g','b','c','m','y','k'];
                            %linew= 1.6;
                            color=0;
							segmentindicesU= [];
                            hold on
                            for j=1:size(Y,1)
                                currclusternumber = clustersU(j);
                                if j==1
                                    oldclusternum = clustersU(j);
                                elseif(oldclusternum ~= currclusternumber)
                                    oldclusternum = currclusternumber;
                                    color=color+1;
									segmentindicesU = [j,segmentindicesU];
                                end
                                %if color>=N
                                 %   linew=linew+.2;
                                %end
								%disp(color);
								coloridx = mod(color,N);
								if coloridx==0
                                    coloridx=N;
                                end
                                %coloridx=colorVar(color);
                                plot(Y(j,:),'color',tempvar(coloridx,:),'Linewidth',1.6); %plot county timeseries.
                                
                            end
                            %{
                            for j=1:size(Y,1)
                                coloridx=mod(j,N);
                                if coloridx==0
                                    coloridx=N;
                                end
                                plot(Y(j,:),'color',tempvar(coloridx,:),'Linewidth',1.6); %plot county timeseries.
                            end
                            %}
                           % xlim([0,size(Y,2)]);


                            %get segment boundaries for V
                            segmentindicesV= [];
                            for j=1:size(X,2)
                                currclusternumber = clustersV(j);
                                if j==1
                                    oldclusternum = clustersV(j);
                                elseif(oldclusternum ~= currclusternumber)
                                    oldclusternum = currclusternumber;
                                    segmentindicesV = [j,segmentindicesV];
                                end
                            end
                            %{
                            oldclusternum=-1;
                            %get segment boundaries for U
                            for j=1:size(Y,1)
                                currclusternumber = clustersU(j);
                                if j==1
                                    oldclusternum = clustersU(j);
                                elseif(oldclusternum ~= currclusternumber)
                                    oldclusternum = currclusternumber;
                                    segmentindicesU = [j,segmentindicesU];
                                end
                            end
                            %}
                            resultfilename=strcat('lambda_1_',mat2str(lambda_1),'_lambda_2_',mat2str(lambda_2),...
                             '_lambda_3_',mat2str(lambda_3),'_numiter_',mat2str(maxIterationbCluster),...
                             '_numcluster_',mat2str(nbCluster));  

                            segmentindicesU=sort(segmentindicesU);
                            segmentindicesV=sort(segmentindicesV); %sort in ascending order.
                            s = sprintf('lam1 %d lam2 %d lam3 %d',lambda_1,lambda_2,lambda_3);
                            disp(s);
                            disp('segment indices V, U');
                            disp(size(segmentindicesV));
							disp(size(segmentindicesU));
                            %{
                            %segment
                            segoutputfileV=strcat('osc_segment_indices_V_',resultfilename);
                            segoutputfileU=strcat('osc_segment_indices_U_',resultfilename);
                            segmentfilenameV=strcat(resultsoutputdir,segoutputfileV,'.csv');
                            segmentfilenameU=strcat(resultsoutputdir,segoutputfileU,'.csv');
                            csvwrite(segmentfilenameV,segmentindicesV);
                            csvwrite(segmentfilenameU,segmentindicesU);
                            
                            %V matrix write
                            
                            vmatrixoutputfile=strcat(resultsoutputdir,'V_matrix_',resultfilename,'.csv');
                            fprintf('vmatrixoutputfile %s\n',vmatrixoutputfile);
                            dlmwrite(vmatrixoutputfile,V);
                            
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

                            %figure output name UV
                            
                            
                            fignameV=strcat('osc_segments_UV_',resultfilename);
                            figurefilenameV=strcat(resultsoutputdir,fignameV,'.pdf');
                
                            y1=get(gca,'ylim');
                            x1=get(gca,'xlim');
                            %clustercolor=['b','y','r','k','g'];
                            %clustercolor=['y','y','y','y','y'];
                            for j=1:size(segmentindicesV,2)
                                line([segmentindicesV(j),segmentindicesV(j)],y1,'Color','k');
                            end
							%}
                            %{
                            alphaval=0.00001;
                            for j=1:size(segmentindicesV,2)-1
                                startidx=segmentindicesV(j);
                                fprintf('Startidx %d and cluster %d\n',startidx,clustersV(startidx));
                                endidx=segmentindicesV(j+1);
                                clusternum=mod(clustersV(startidx),size(clustercolor,2));
                                if clusternum==0
                                    clusternum=size(clustercolor,2);
                                end
                                color=clustercolor(clusternum);
                                p=patch([startidx endidx endidx startidx],[1 1 y1(2) y1(2)],color);
                                p.FaceAlpha=alphaval;
                            end
                            
                            %populate the area before the first segment
                            endidx=segmentindicesV(1); %get the end idx of the first segment.
                            xlimit=x1(1); %get the beginning of the xlimit of the plot
                            clusternum=clustersV(1); %get the cluster number of the first timestep.
                            color=clustercolor(clusternum);
                            p=patch([xlimit,endidx,endidx,xlimit],[1,1,y1(2),y1(2)],color);
                            p.FaceAlpha=alphaval;
                            
                            %populate the end segment
                            endidx=segmentindicesV(end); %get the end idx of the last segment.
                            xlimit=x1(2);
                            clusternum=clustersV(endidx);
                            color=clustercolor(clusternum);
                            p=patch([endidx,xlimit,xlimit,endidx],[1,1,y1(2),y1(2)],color);
                            p.FaceAlpha=alphaval;
                            
                            set([gca],'FontSize', 18);
                            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                            saveas(gcf,figurefilenameV);
                            print(figurefilenameV,'-dpdf','-bestfit');
                            hold off
                            %}
                            %{
                            %figure output name U
                            fignameU=strcat('osc_segments_U_',resultfilename);
                            figurefilenameU=strcat(resultsoutputdir,fignameU,'.pdf');

                            segmentindicesU=sort(segmentindicesU); %sort in ascending order.
                            y1=get(gca,'ylim');
                            x1=get(gca,'xlim');
                            clustercolor=['b','y','r','k','g'];
                            alphaval=0.2;
                            for j=1:size(segmentindicesU,2)-1
                                startidx=segmentindicesU(j);
                                fprintf('Startidx %d and cluster %d\n',startidx,clustersU(startidx));
                                endidx=segmentindicesU(j+1);
                                clusternum=mod(clustersU(startidx),size(clustercolor,2));
                                if clusternum==0
                                    clusternum=size(clustercolor,2);
                                end
                                color=clustercolor(clusternum);
                                p=patch([startidx endidx endidx startidx],[1 1 y1(2) y1(2)],color);
                                p.FaceAlpha=alphaval;
                            end

                            %populate the area before the first segment
                            endidx=segmentindicesU(1); %get the end idx of the first segment.
                            xlimit=x1(1); %get the beginning of the xlimit of the plot
                            clusternum=clustersU(1); %get the cluster number of the first timestep.
                            color=clustercolor(clusternum);
                            p=patch([xlimit,endidx,endidx,xlimit],[1,1,y1(2),y1(2)],color);
                            p.FaceAlpha=alphaval;

                            %populate the end segment
                            endidx=segmentindicesU(end); %get the end idx of the last segment.
                            xlimit=x1(2);
                            clusternum=clustersU(endidx);
                            color=clustercolor(clusternum);
                            p=patch([endidx,xlimit,xlimit,endidx],[1,1,y1(2),y1(2)],color);
                            p.FaceAlpha=alphaval;

                            set([gca],'FontSize', 18);
                            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                            saveas(gcf,figurefilenameU);
                            print(figurefilenameU,'-dpdf','-bestfit');
                            hold off
                            %}
                            %{
                            clusterfilenameU=strcat(resultsoutputdir,'clusters_U',resultfilename,'.csv');
                            fileID = fopen(clusterfilenameU,'w');
                            fprintf(fileID,'%d\n',clustersU);
                            fclose(fileID);
                            clusterfilenameV=strcat(resultsoutputdir,'clusters_V',resultfilename,'.csv');
                            fileID = fopen(clusterfilenameV,'w');
                            fprintf(fileID,'%d\n',clustersV);
                            fclose(fileID);
                            %}
                            close(gcf); 
                    end
                end
		    end
        end
    end
	rmpath(paths);
	rmpath(paths2);
end
