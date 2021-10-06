%  Hurricane demo
%If you want to convert this into a function, the following variables need
%to be set
%inputfile, visfile,maxIterationbCluster,lambda_1,lambda_2,gamma_1,gamma_2
function [] = demo_chickendance_seg_exog3_fast(resultsoutputdir,inputfile,laplacefile)
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
        %chickendance
        lam1=[-10,-12];
        lam2=[1.5,2,2.5,3,3.5];
        lam3 =[1.5,2,2.5,3,3.5];
		beta=[10];
		%}
       
        lam1=[-1];%-10
        lam2=[2];
        lam3=[3];%3.5
        beta=[1];%10
        
        %nbclustersV=[2];
        %maxiter=[300];
    iter = 300;
	%for nbClusterV=nbclustersV
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
                        %[V, U, funVal,it] = OSC_Exog3(X, L, 2, lambda_1, lambda_2, lambda_3, gamma_1, gamma_2, gamma_3, p,... 
                         %                    beta_1, maxIterationbCluster, diag, diag);
						[V, R, funVal, it] = OSC_fast(X, lambda_1, lambda_2, gamma_1, gamma_2, p, maxIterationbCluster, diag);
						fprintf('param:%.2f,%.2f,%.2f,%.2f,%.3f converged %d\n',lambda_1,lambda_2,lambda_3,...
                                    beta_1,funVal(it),it);
						resultfilename=strcat('lam1_',mat2str(lambda_1),'_lam2_',mat2str(lambda_2),'_lam3_',...
                                mat2str(lambda_3),'_beta_',mat2str(beta_1)); 
						%V matrix write
                        vmatrixoutputfile=strcat(resultsoutputdir,'V_matrix_',resultfilename,'.csv');
                        fprintf('vmatrixoutputfile %s\n',vmatrixoutputfile);
                        dlmwrite(vmatrixoutputfile,V);
                                 
                        %U matrix write
                        %umatrixoutputfile=strcat(resultsoutputdir,'U_matrix_',resultfilename,'.csv');
                        %fprintf('umatrixoutputfile %s\n',umatrixoutputfile);
                        %dlmwrite(umatrixoutputfile,U);
					end
				end
			end
		end
	rmpath(paths);
	rmpath(paths2);
end
						
