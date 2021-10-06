function [ V, R, funVal, iteration ] = OSCExog_fast( X,L, lambda_1, lambda_2, gamma_1, gamma_2, p, maxIterations, diagconstraint)
%% Problem
% Modify this to add exogenous inputs.
%
% min L(Z) = 1/2 ||X - UXV||^2_F + lambda_1 * ||U||_1 + (beta_1)/2 Tr(U'LU) 
%                + lambda_2 * ||V||_1 + lambda_3 * ||VR||_{1,2} 
% 
% where ||B||_2/1 = ||b_1||_2 + ||b_2||_2 + ... + ||b_n||_2
%
%% Solution 
% We solve this problem via the ADMM (Alternating direction method of multipliers) variant
% of Augmented Lagrangian method as follows:
% 
% Let K = V and J = U , P = KR then we have
% 
%
% T(U,V,J,K,P) = 1/2||X - JXK||^2_F + lambda_1 * ||U||_1 + (beta_1)/2 * Tr(J'LJ) 
%              + lambda_2 * ||V||_1  +lambda_3 * ||P||_{1,2}  
%              + <G,V - K> + (gamma_1)/2 * ||V - K||^2_F 
%              + <H,U - J> + (gamma_2)/2 * ||U - J||^2_F
%              + <F,P - KR>+ (gamma_3)/2 * ||P - KR||_F^2
%

%% Matrix Dimensions

% n = Number of Counties
% m = Number of Timesteps
% dim(X) = n X m

% dim(J) = n X n
% dim(U) = n X n
% dim(H) = n X n
% dim(L) = n X n This is the Laplacian matrix.

% dim(V) = m X m
% dim(G) = m X m
% dim(K) = m X m

% dim(R) = m X m - 1
% dim(P) = m X m - 1
% dim(F) = m X m - 1

%% 
if (~exist('diagconstraint','var'))
    diagconstraint = 0;
end

funVal = zeros(maxIterations,1);
%rename xn to n and m = number of counties.
[n, m, ~] = size(X);  % Num Counties X Num Timesteps

%Temporal Matrices
V = zeros(m, m);
G = zeros(m, m);
K = zeros(m, m); % K = V

R = (triu(ones(m,m-1),1) - triu(ones(m, m-1))) + (triu(ones(m, m-1),-1)-triu(ones(m, m-1)));
R = sparse(R);
P = zeros(m, m-1);
F = zeros(m, m-1); 

%Spatial Matrices
U = zeros(n,n);
J = zeros(n,n);
H = zeros(n,n);

for iteration=1:maxIterations
    fprintf('Iteration Number %d \n',iteration);
    
    %% Step 1  V Update. 
    DV = K - (G/gamma_1);
    [vm, vn, ~] = size(DV);
    rolled_v = reshape(DV,vm*vn,1);

    rolled_z = shrink_l1(rolled_v, lambda_2/gamma_1);

    V = reshape(rolled_z, vm, vn);

    % Set V diag to 0
    if (diagconstraint)
        V(logical(eye(size(V)))) = 0;
    end

    %% Step 2 K Update.
    A = X'*(J'*J)*X + gamma_1*speye(m,m);
    B = gamma_3*(R*R');
    C = -(X'*J*X + G + gamma_1*V + gamma_3*P*R'  + F*R');
    K = lyap(A, B, C);

    %% Step 3 P Update.
    DV = K*R - (1/gamma_3)*F;
    P = mysolve_l1l2(DV, lambda_3/gamma_3);

    %% Step 4 U Update.
    DU = J - (H/gamma_2);
    [um, un, ~] = size(DU);
    rolled_u = reshape(DU,um*un,1);
    rolled_z = shrink_l1(rolled_u, lambda_1/gamma_2);
    U = reshape(rolled_z, um, un);

    % Set U diag to 0
    if (diagconstraint)
        U(logical(eye(size(U)))) = 0;
    end
    
    %% Step 5 J Update.
    A = beta_1*L*J;
    B = X*(K*K')*X + gamma_2*speye(n,n);
    C = X*K*X' + H + gamma_2*U;
    J = lyap(A,B,C);
    
    %% Step 6

    G = G + gamma_1 * (V - K);

    %% Step 7

    F = F + gamma_2 * (P - K*R);

    %% Step 8

    gamma_1 = p * gamma_1;
    gamma_2 = p * gamma_2;

    %% Calculate function values
    funVal(iteration) = .5 * norm(X - X*V,'fro')^2 + lambda_1*norm(V,1) +lambda_2*l2l1norm(V*R);

    if iteration > 1
        if funVal(iteration) < 1*10^-3
            break
        end
    end

    if iteration > 100
        if funVal(iteration) < 1*10^-3 || funVal(iteration-1) == funVal(iteration) ...
                || funVal(iteration-1) - funVal(iteration) < 1*10^-3
            break
        end
    end


end

end