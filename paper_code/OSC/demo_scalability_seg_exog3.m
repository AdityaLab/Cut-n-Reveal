%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_scalability_seg_exog3(resultsoutputdir,inputfile,laplacefile)
    paths = genpath('libs/ncut');
	paths2 = genpath('libs/linspecer');
    fprintf('results dir %s\n',resultsoutputdir);
	addpath(paths);
	addpath(paths2);

	%% Load Hurricane Data

	rng(1);

	X = csvread(inputfile,1);
	Y = csvread(inputfile,1); 

	X=X'; %Size(X) will be equal to 366 X 625
	Y=Y'; %This is the data to visualize 366 X 625.
	corruption = 0;

	w = randn(size(X)) * corruption;
	X = X + w;

	%Hyperparameter.
		%{
        %NILM
        lam1=[-10,0.7,10];%-15,-5
        lam2=[-10,0.2,0.7,10];%5,10,15
        lam3=[8,10,12];%0.1
        %}
        %beta=[10];
        %beta=[2,10,20,50];
        %nbclustersV=[2,3,4];
       
        
        lam1=[-10];
        lam2=[10];
        lam3=[10];
        beta=[1];
        
        nbclustersV=[3];
        %maxiter=[300];
    iter = 300;
    %lossfile= strcat(resultsoutputdir,'loss4.txt');
	%lossID = fopen(lossfile,'w');
    %fprintf(lossID,'lambda1, lambda2, lambda3, beta1, gam1, gam2, gam3, clusterV, clusterU, loss\n');
    timefile = strcat(resultsoutputdir,'nilm_500timesteps.csv');
	%timefile = strcat(resultsoutputdir,'nilm_480timeseries.csv.csv');
	time = fopen(timefile,'w');
    ttl_time=0;
	fprintf('start\n');
	for nbClusterV=nbclustersV
		%for nbClusterU=nbclustersU
        for lambda_1=lam1 
			for lambda_2=lam2
                for lambda_3=lam3
                    for beta_1=beta
                        %% OSC
                        maxIterationbCluster = iter;
                        gamma_1 = 0.01;
                        gamma_2 = 0.01;
                        gamma_3 = 0.01;
                        p = 1.1;
                        diag = 1;
                        L = importdata(laplacefile);
                        %k=[5,10];
                        tic;
                        [V, U, funVal,it] = OSC_Exog3(X, L, 2, lambda_1, lambda_2, lambda_3, gamma_1, gamma_2, gamma_3, p,... 
                                             beta_1, maxIterationbCluster, diag, diag);
                        et=toc;
						fprintf(time,'%.4f\n',et);
                        ttl_time+=time;
						
						vAffine = V'*V;
                        vInval = nnz(isnan(vAffine))+nnz(isinf(vAffine))+all(vAffine(:)==0);
						if funVal(it)>=0 && vInval==0
                            fprintf('Iteration Number %d \n',it);
                            [oscVclusters,oscVeigenvectors,oscVeigenvalues] = ncutW(vAffine,nbClusterV); %normalized V cuts.
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
                            %{
                            s = sprintf('lam1 %.3f lam2 %.3f lam3 %.3f beta %.2f %.3f',lambda_1,lambda_2,lambda_3,beta_1,funVal(it));
                            disp(s);
                            disp('segment indices V');
                            disp(size(segmentindicesV));
							%disp(size(segmentindicesU));
							fprintf(lossID,'%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%.3f\n',lambda_1,lambda_2,lambda_3,...
                                    beta_1,gamma_1,gamma_2,gamma_3,nbClusterV,funVal(it));
							fV = prune_candidate_segment(.05,segmentindicesV,nbClusterV,size(X,2),0); 
                            %fU = prune_candidate_segment(.05,segmentindicesU,nbClusterU,size(X,1),0);
                            %fV=1;8,25
                            if fV & size(segmentindicesV,2)>=0 & size(segmentindicesV,2)<=10 
                                resultfilename=strcat('lam1_',mat2str(lambda_1),'_lam2_',mat2str(lambda_2),'_lam3_',...
                                mat2str(lambda_3),'_clusV_',mat2str(nbClusterV),'_beta_',mat2str(beta_1));
                                %clustering of non-normalized timeseries.
                                %segment
                                 segoutputfileV=strcat('segV_',resultfilename);
                                 %segoutputfileU=strcat('segU_',resultfilename);
                                 segmentfilenameV=strcat(resultsoutputdir,segoutputfileV,'.csv');
                                 %segmentfilenameU=strcat(resultsoutputdir,segoutputfileU,'.csv');
                                 csvwrite(segmentfilenameV,segmentindicesV);
                                 %csvwrite(segmentfilenameU,segmentindicesU);
                                 
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
                                 dlmwrite(affinitymatrixoutputfileV,(V'*V));
                                 %{
                                 %Affinity matrix U write
                                 affinitymatrixoutputfileU=strcat(resultsoutputdir,'U_affinity_matrix_',resultfilename,'.csv');
                                 dlmwrite(affinitymatrixoutputfileU,(U*U'));
                                 %}
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
                            %}
                        end
						
                   end
                end
		    end
        end
    end
    %fclose(lossID);
	fclose(time);
	rmpath(paths);
	rmpath(paths2);
end
