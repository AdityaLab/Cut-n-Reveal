function [ E ] = find_exp(alpha, lmda, save_dir)
%FIND_EXP Summary of this function goes here
%   given a segmentation of the segmentation of the data (the difference
%   between the adjacency time segments), the laplacian matrix among the
%   time series, the weights for the regularizations (alpha for the
%   laplacian regularization, lmda for the lasso regularization); find the
%   explanation for all the cut points. 
%   The network where the laplacian matrix comes from should be a connected
%   network.
    %Below is a toy dataset just for testing.
    %diff = [1,3,2,1,4;3,3,1,1,1;1,0,5,4,0];
    %L = [2,-1,-1,0,0;
    %	-1,1,0,0,0;
    %	-1,0,2,-1,0;
    %	0,0,-1,2,-1;
    %	0,0,0,-1,1];
    %alpha = 0.2;
    %lmda = 0.2;
    load(strcat(save_dir, 'dif.mat'));
    load(strcat(save_dir, 'L.mat'));
    dif = double(dif);
    L = double(L);
    H = alpha * L * 2;
    lb = zeros(size(dif,2), 1);
    up =  ones(size(dif,2), 1);
    Aeq = ones(1, size(dif,2));
    beq = [1];
    E = [];
    for ii = 1:size(dif, 1)
        f = (-1 * (dif(ii,:) - lmda))';
        [x, fval, exitflag, output, lambda] = quadprog(H, f, [], [], Aeq, beq, lb, up);
        E = [E;x'];
    end
    E = E';
    fid = fopen(strcat(save_dir, 'E_',mat2str(alpha),'.txt'), 'wt');
    for ii = 1: size(E, 1)
        fprintf(fid, '%20.18f\t', E(ii,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
    %dif = double(dif);
    %L = double(L);
    %J = ones(size(dif,1), size(dif,2)).*lmda;
    %invL = pinv(L);
    %E = (invL * (dif.' - J.')) ./alpha;
end

