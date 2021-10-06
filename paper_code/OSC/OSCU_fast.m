function[V, U, funVal, iteration] = OSCU_fast(X, L, lambda_1, lambda_2, lambda_3, gamma_1, gamma_2,... 
                                                gamma_3, beta_1, p, maxIterations, diag1, diag2)
%% Problem
% X n*m matrix
% min U,V = 1/2 ||X - UXV||^2_F + lambda_1 * ||U||_1 + gamma1/2 * Tr(U'LU) +
%           + lambda_2 * ||V||_1 + lambda_3 * ||ZR||_2/1
% 
% where ||B||_2/1 = ||b_1||_2 + ||b_2||_2 + ... + ||b_n||_2
%
%% Solution 
% We solve this problem via the ADMM (Alternating direction method of multipliers) variant
% of Augmented Lagrangian method as follows:
% 
% Let J = U, K =V, P = KR then we have
% 
% L(U,V,J,K,P) = 1/2||X - JXK||^2_F + lambda_1 * ||U||_1 + gamma1/2 * Tr(J'LJ)
%                lambda_2 * ||V||_1 + lambda_3 * ||P||_2/1 + <G,V - K>
%               + gamma_1/2 * ||V - K||^2_F + <H,U - J> + gamma_2/2 * ||U - J||^2_F 
%               + <F,P - KR> + gamma_3/2 * ||P - KR||^2_F 
% 

if (~exist('diag1','var'))
    diag1 = 0;
end

if (~exist('diag2','var'))
    diag2 = 0;
end
funVal = zeros(maxIterations,1);
%Initialization
[~, xn, ~] = size(X);
[~, ln, ~] = size(L);
U = zeros(ln,ln);
V = zeros(xn,xn);

R = (triu(ones(xn,xn-1),1) - triu(ones(xn, xn-1))) + (triu(ones(xn, xn-1),-1)-triu(ones(xn, xn-1)));
R = sparse(R);

J = zeros(ln,ln); %J=U
K = zeros(xn,xn); %K=V
P = zeros(xn, xn-1); %P=KR

G = zeros(xn,xn);
H = zeros(ln,ln);
F = zeros(xn,xn-1);

for iteration=1:maxIterations
    %step 1 solve for V
    Z = K - (G/gamma_1);
    [zm, zn, ~] = size(Z);
    rolled_z = reshape(Z,zm*zn,1);

    rolled_v = shrink_l1(rolled_z, lambda_2/gamma_1);

    V = reshape(rolled_v, zm, zn);

    % Set Z diag to 0
    if (diag1)
        V(logical(eye(size(V)))) = 0;
    end
    %step2 solve for U
    Y = J - (H/gamma_2);
    [ym, yn, ~] = size(Y);
    rolled_y = reshape(Y,ym*yn,1);
    rolled_u = shrink_l1(rolled_y, lambda_1/gamma_2);
    U = reshape(rolled_u,ym,yn);
    if (diag2)
        U(logical(eye(size(U)))) = 0;
    end
    
    %step 3 solve for K
    A_1 = X'*J'*J*X + gamma_1 * speye(xn,xn);
    B_1 = gamma_3 * R * R';
    C_1 = -(X'*J*X + G + lambda_1*V + F*R' + lambda_3 * P * R');
    K = lyap(A_1,B_1,C_1);
    
    %step 4 solve for J
    A_2 = beta_1 * L;
    B_2 = X * K * K' * X' + gamma_2 * speye(ln,ln);
    C_2 = X * K * X' + H + gamma_2 * U;
    J =lyap(A_2,B_2,C_2);
    
    %step 5 solve for P
    M = K*R - (1/gamma_3)*F;
    P = mysolve_l1l2(M, lambda_3/gamma_3);
    
    %step 6
    G = G + gamma_1 * (V-K);
    %step 7
    H = H + gamma_2 * (U-J);
    %step 8
    F = F + gamma_3 * (P-K*R);
    %step 9
    gamma_1 = p*gamma_1;
    gamma_2 = p*gamma_2;
    gamma_3 = p*gamma_3;
    
     %% Calculate function values
    funVal(iteration) = .5 * norm(X - U*X*V,'fro')^2 + lambda_1*norm(U,1) +...
                         beta_1/2 * trace(U'* L * U) + lambda_2*norm(V,1)+lambda_3*l2l1norm(V*R); 

    if iteration > 1
        if funVal(iteration) < 1*10^-3
            break
        end
    end
    %fprintf('\nfunction value %.12f, iteration %d\n',funVal(iteration),iteration);
    if iteration > 100
        if funVal(iteration) < 1*10^-4 || funVal(iteration-1) == funVal(iteration) ...
                || funVal(iteration-1) - funVal(iteration) < 1*10^-4
            break
        end
    end
end


end
