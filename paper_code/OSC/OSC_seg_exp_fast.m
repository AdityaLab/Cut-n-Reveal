function [ V, R, funVal, iteration ] = OSC_seg_exp_fast(X,S, lambda_1, beta_2, beta_1, p, maxIterations, diagconstraint)
%% Problem
%
% min L(V) = 1/2 ||X - XZ||^2_F + lambda_1 * ||Z||_1 + \beta_1/2||VR -
% S||^2_F  
%  S = m X m - 1 matrix.
%  X = n X m  data matrix
%  
%
%
%% Solution 
% We solve this problem via the ADMM (Alternating direction method of multipliers) variant
% of Augmented Lagrangian method as follows:
% 
% Let V = K and let S denote the explanation matrix.
% 
% T(V,K) = 1/2||X - XV||^2_F + lambda_1 * ||V||_1
% + beta_1/2 * ||KR - S||^2_F + <G,V - K> + beta_2/2 * ||V - K||^2_F
%

if (~exist('diagconstraint','var'))
    diagconstraint = 0;
end

funVal = zeros(maxIterations,1);

[~, m, ~] = size(X); %xn is number of time steps

K = zeros(m, m); % K = V
R = (triu(ones(m,m-1),1) - triu(ones(m, m-1))) + (triu(ones(m, m-1),-1)-triu(ones(m, m-1)));
R = sparse(R);

G = zeros(m, m);

V = zeros(m, m);

for iteration=1:maxIterations
    tic
    fprintf('Iteration Number %d \n',iteration);
    %% Step 1 V update. Need to check whether this is correct or not.
    V = K - (G/beta_2);
    [vm, vn, ~] = size(V);
    rolled_v = reshape(V,vm*vn,1);

    rolled_z = shrink_l1(rolled_v, lambda_1/beta_2);

    V = reshape(rolled_z, vm, vn);

    % Set Z diag to 0
    if (diagconstraint)
        V(logical(eye(size(V)))) = 0;
    end
    
    %% Step 2  K update.
    A = X'*X + beta_2*speye(m,m);
    B = beta_1*(R*R');
    C = -(X'*X + beta_1*(S*R') + G + beta_2*V);

    K = lyap(A, B, C); %this is the step that takes the most time.
    
    %% Step 3
    
    G = G + beta_2 * (V - K);
    
    %% Step 5
    
    beta_2 = p * beta_2;
    beta_1 = p * beta_1;
    %% Calculate function values
    
    funVal(iteration) = .5 * norm(X - X*V,'fro')^2 + lambda_1*norm(V,1) +(beta_1/2)*norm(V*R - S,'fro')^2;
   
    if iteration > 1
        if funVal(iteration) < 1*10^-3
            break
        end
    end
    
    if iteration > 100
        if funVal(iteration) < 1*10^-4 || funVal(iteration-1) == funVal(iteration) ...
                || funVal(iteration-1) - funVal(iteration) < 1*10^-3
            break
        end
    end

toc
end

end