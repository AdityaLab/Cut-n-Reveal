%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_hurricane_seg_exog3_fast(resultsoutputdir,inputfile,visfile,laplacefile)
    %resultsoutputdir = '../result/Harvey_exog/';
    %inputfile = '../data/Harvey_60min_sample_normalized.csv';
    %visfile = '../data/Harvey_60min_sample.csv';
    paths = genpath('libs/ncut');
	paths2 = genpath('libs/linspecer');
    fprintf('results dir %s\n',resultsoutputdir);
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
       %%{
        lam1=[0.1,0.5,0.7,2];
        lam2=[0.1,0.5,0.7,2];
        lam3 =[0.1,0.5,0.7,2];
        beta=[0.7,1];
        %%}
        %beta=[2,10,20,50];
        %nbclustersU=[2,3,4];
        %nbclustersV=[2,3,4];
        %{
        lam1=[0.1,2];
        lam2=[0.1,2];
        lam3=[0.1,2];
        beta=[0.7,1];
		k=5;
        %}
        %nbclustersU=[2];
        %nbclustersV=[2,3,4];
        %maxiter=[300];
    iter = 1000;
	%lam1=[0.1];
    %lam2=[0.1];
    %lam3=[0.1,0.3,0.5];
    %beta=[0.7,1];
	nbclustersV=[3];
    %hurricane Matthew params
    %{lam1=[0.1];
    lam2=[0.1];
    lam3=[0.1];
    beta=[1];
    nbclustersV=[4];
    %}
	k=5;
    lossfileName = strcat(resultsoutputdir,'loss.txt');
	lossID = fopen(lossfileName,'w');
    fprintf(lossID,'lambda1, lambda2, lambda3, beta1, gam1, gam2, gam3, clusterV, clusterU, loss\n');
    for nbClusterV=nbclustersV
		%for nbClusterU=nbclustersU
            for lambda_1=lam1 
			    for lambda_2=lam2
                    for lambda_3=lam3
                        for beta_1=beta
                            %% OSC
                            maxIterationbCluster = iter;
                            gamma_1 = 0.5;
                            gamma_2 = 0.5;
                            gamma_3 = 0.5;
                            p = 1.1;
                            diag = 1;
                            L = importdata(laplacefile);
                            %k=[5,10];
                            tic;
                            [V, U, funVal,it] = OSC_Exog3(X, L, k, lambda_1, lambda_2, lambda_3, gamma_1, gamma_2, gamma_3, p,... 
                                             beta_1, maxIterationbCluster, diag, diag);
                            toc;
                            fprintf('Iteration Number %d \n',it);
                            [oscVclusters,oscVeigenvectors,oscVeigenvalues] = ncutW(V'*V,nbClusterV); %normalized V cuts.
                            %[oscUclusters,oscUeigenvectors,oscUeigenvalues] = ncutW(U*U',nbClusterU); %normalized U cuts
                            clustersV = denseSeg(oscVclusters,1);
                            %clustersU = denseSeg(oscUclusters,1);
                            %disp('V cluster matrix');
                            %disp(clustersV);                           
                            %disp('U cluster matrix');
                            %disp(clustersU);
                            %segmentindicesU= [];
                            segmentindicesV= [];
                            oldclusternum=-1;
                            for j=1:size(X,2)
                                currclusternumber = clustersV(j);
                                if j==1
                                    oldclusternum = clustersV(j);
                                elseif(oldclusternum ~= currclusternumber)
                                    oldclusternum = currclusternumber;
                                    segmentindicesV = [j,segmentindicesV];
                                end
                            end
                            segmentindicesV=sort(segmentindicesV); %sort in ascending order.
                            s = sprintf('lam1 %.3f lam2 %.3f lam3 %.3f beta %.2f %.3f',lambda_1,lambda_2,lambda_3,beta_1,funVal(it));
                            disp(s);
                            disp('segment indices V');
                            disp(size(segmentindicesV));
							%disp(size(segmentindicesU));
                            fprintf(lossID,'%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%.3f\n',lambda_1,lambda_2,lambda_3,...
                                    beta_1,gamma_1,gamma_2,gamma_3,nbClusterV,funVal(it));
                            fV = prune_candidate_segment(.05,segmentindicesV,nbClusterV,size(X,2),0); 
                            %fU = prune_candidate_segment(.05,segmentindicesU,nbClusterU,size(X,1),0);
                            %fV=1;
                            if fV & size(segmentindicesV,2)<10
                                resultfilename=strcat('lam1_',mat2str(lambda_1),'_lam2_',mat2str(lambda_2),'_lam3_',...
                                mat2str(lambda_3),'_clusV_',mat2str(nbClusterV),'_l_',mat2str(k));%,'_clusU_',mat2str(nbClusterU));  
                                %clustering of non-normalized timeseries.
                                %segment
                                 segoutputfileV=strcat('segV_',resultfilename);
                                 %segoutputfileU=strcat('segU_',resultfilename);
                                 segmentfilenameV=strcat(resultsoutputdir,segoutputfileV,'.csv');
                                 %segmentfilenameU=strcat(resultsoutputdir,segoutputfileU,'.csv');
                                 csvwrite(segmentfilenameV,segmentindicesV);
                                 %csvwrite(segmentfilenameU,segmentindicesU);
                                 %{
                                 %V matrix write
                                 vmatrixoutputfile=strcat(resultsoutputdir,'V_matrix_',resultfilename,'.csv');
                                 fprintf('vmatrixoutputfile %s\n',vmatrixoutputfile);
                                 dlmwrite(vmatrixoutputfile,V);
                                 
                                 %U matrix write
                                 umatrixoutputfile=strcat(resultsoutputdir,'U_matrix_',resultfilename,'.csv');
                                 fprintf('umatrixoutputfile %s\n',umatrixoutputfile);
                                 dlmwrite(umatrixoutputfile,U);
                                 %}
                                 %Affinity matrix V write
                                 affinitymatrixoutputfileV=strcat(resultsoutputdir,'V_affinity_matrix_',resultfilename,'.csv');
                                 dlmwrite(affinitymatrixoutputfileV,(V'*V));
                                 
                                 %Affinity matrix U write
                                 affinitymatrixoutputfileU=strcat(resultsoutputdir,'U_affinity_matrix_',resultfilename,'.csv');
                                 dlmwrite(affinitymatrixoutputfileU,(U*U'));
                                 
                                 %figure output name UV      
                                 fignameV=strcat('oscV_',resultfilename);
                                 figurefilenameV=strcat(resultsoutputdir,fignameV,'.pdf');
                                 figure;
                                 N=250;
                                 tempvar = linspecer(N,'qualitative'); %get C for up to 12 different colors.
                                 %colorVar=['r','g','b','c','m','y','k'];
                                 hold on
                                 for j=1:size(Y,1)
                                     coloridx=mod(j,N);
                                     if coloridx==0
                                         coloridx=N;
                                     end
                                     plot(Y(j,:),'color',tempvar(coloridx,:),'Linewidth',1.3); %plot county timeseries.
                                 end
                                 xlim([0,size(Y,2)]); 
                                
                                 y1=get(gca,'ylim');
                                 x1=get(gca,'xlim');
                                 %clustercolor=['b','y','r','k','g'];
                                 %clustercolor=['y','y','y','y','y'];
                                 for j=1:size(segmentindicesV,2)
                                     line([segmentindicesV(j),segmentindicesV(j)],y1,'Color','k');
                                 end    
                                 set([gca],'FontSize', 18);
                                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                                 saveas(gcf,figurefilenameV);
                                 print(figurefilenameV,'-dpdf');
                                 hold off
                                 %{
                                 clusterfilenameU=strcat(resultsoutputdir,'clusters_U_',resultfilename,'.csv');
                                 fileID = fopen(clusterfilenameU,'w');
                                 fprintf(fileID,'%d\n',clustersU);
                                 fclose(fileID);
                                 
                                 clusterfilenameV=strcat(resultsoutputdir,'clusters_V',resultfilename,'.csv');
                                 fileID = fopen(clusterfilenameV,'w');
                                 fprintf(fileID,'%d\n',clustersV);
                                 fclose(fileID);
                                 %}
                                 close(gcf);
                                 fprintf('plotting figure\n');
                            end
                        end
                   end
                end
		    end
        %end
    end
    fclose(lossID);
	rmpath(paths);
	rmpath(paths2);
end
