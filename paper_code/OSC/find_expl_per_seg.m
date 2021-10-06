function segmentindicesV = find_expl_per_seg(resultsoutputdir,X,laplacefile,segNo)
	paths = genpath('libs/ncut');
	addpath(paths);
	%Hyperparameter
     %{
     lam1=[0.1,0.5,0.7,2];
     lam2=[0.1,0.5,0.7,2];
     lam3 =[0.1,0.5,0.7,2];
     beta=[0.7,1];   
     %}
     lam1=[0.1];
     lam2=[0.02];
     lam3=[0.02];
     beta=[0.5];
     
    %nbClusterU=3;
    nbClusterV=2;
    iter = 300;
    %lossfileName = strcat(resultsoutputdir,'loss.txt');
	%lossID = fopen(lossfileName,'w');
    %fprintf(lossID,'lambda1, lambda2, lambda3, beta1, gam1, gam2, gam3, clusterV, clusterU, loss\n');
    for lambda_1=lam1 
		for lambda_2=lam2
            for lambda_3=lam3
                for beta_1=beta
                    %OSC
                    maxIterationbCluster = iter;
                    gamma_1 = 0.1;
                    gamma_2 = 0.1;
                    gamma_3 = 0.1;
                    p = 1.1;
                    diag = 1;
                    L = importdata(laplacefile);
                    [V, U, funVal,it] = OSC_Exog3(X, L, 3, lambda_1, lambda_2, lambda_3, gamma_1, gamma_2, gamma_3, p,... 
                                             beta_1, maxIterationbCluster, diag, diag);
                    fprintf('Iteration Number %d \n',it);
                    [oscVclusters,oscVeigenvectors,oscVeigenvalues] = ncutW(V'*V,nbClusterV); %normalized V cuts
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
                    s = sprintf('lam1 %.3f lam2 %.3f lam3 %.3f beta %.2f segSz %d segNo %d\n',lambda_1,lambda_2,lambda_3,beta_1,size(segmentindicesV),segNo);
                    %disp(s);
                    %fprintf(lossID,'%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%.3f\n',lambda_1,lambda_2,lambda_3,...
                    %beta_1,gamma_1,gamma_2,gamma_3,nbClusterV,funVal(it));
                    fV = prune_candidate_segment(.05,segmentindicesV,nbClusterV,size(X,2),0); 
                    %fV=1;
                    if fV && size(segmentindicesV,2)<3   
                                 %{
                                 %V matrix write
                                 vmatrixoutputfile=strcat(resultsoutputdir,'V_matrix_',resultfilename,'.csv');
                                 fprintf('vmatrixoutputfile %s\n',vmatrixoutputfile);
                                 dlmwrite(vmatrixoutputfile,V);
                                 
                                 %U matrix write
                                 umatrixoutputfile=strcat(resultsoutputdir,'U_matrix_',str(segNo),'.csv');
                                 fprintf('umatrixoutputfile %s\n',umatrixoutputfile);
                                 dlmwrite(umatrixoutputfile,U);
                                 
                                 %Affinity matrix V write
                                 affinitymatrixoutputfileV=strcat(resultsoutputdir,'V_affinity_matrix_',resultfilename,'.csv');
                                 dlmwrite(affinitymatrixoutputfileV,(V'*V));
                                 %}
                                 %Affinity matrix U write
                                 affinitymatrixoutputfileU=strcat(resultsoutputdir,'U_affinity_matrix_',mat2str(segNo),'.csv');
                                 dlmwrite(affinitymatrixoutputfileU,(U*U'));
                                %{ 
                                 clusterfilenameV=strcat(resultsoutputdir,'clusters_V',resultfilename,'.csv');
                                 fileID = fopen(clusterfilenameV,'w');
                                 fprintf(fileID,'%d\n',clustersV);
                                 fclose(fileID);
                                 %}
                                 return;
                    end
                end
            end
        end
	end
    %fclose(lossID);
    segmentindicesV=[];
	rmpath(paths);
end