function [ U_idx, funValU, iteration ] = OSC_Exog2( Xu, L, lambda_1, gamma_1, beta_1, p, maxIterations, diag)
%% Problem
% X n*m matrix
% % min U= 1/2 ||X - UX||^2_F + lambda_1 * ||U||_1 + beta1/2 * Tr(U'LU)
%% Solution 
% We solve this problem via the ADMM (Alternating direction method of multipliers) variant
% of Augmented Lagrangian method as follows:
% 
% Let J = U
% L(U,J) = 1/2||X - JX||^2_F + lambda_1 * ||U||_1 + beta_1/2 * Tr(J'LJ)
%                + <G,U - J> + gamma_1/2 * ||U - J||^2_F 
%
if (~exist('diag','var'))
    diag = 0;
end 
funValU = zeros(maxIterations,1);
%[~, xn, ~] = size(Xu);
[~, ln, ~] = size(L);
U_idx = zeros(ln,ln);


J = zeros(ln,ln); %J=U
G = zeros(ln,ln); %lagrange

for iteration=1:maxIterations
    %step 1 solve solve for U
    Z = J - (G/gamma_1);
    [zm, zn, ~] = size(Z);
    rolled_z = reshape(Z,zm*zn,1);

    rolled_u = shrink_l1(rolled_z, lambda_1/gamma_1);

    U_idx = reshape(rolled_u, zm, zn);

    % Set Z diag to 0
    if (diag)
        U_idx(logical(eye(size(U_idx)))) = 0;
    end
    %step 2 solve for J
    A = beta_1 * L;
    B = Xu*Xu' + gamma_1 *speye(ln,ln);
    C = -(Xu*Xu' + G + gamma_1*U_idx);
    J = lyap(A,B,C);
    
    %step 3
    G = G + gamma_1 * (U_idx-J);
    %step 4
    gamma_1 = p*gamma_1;
    
     %% Calculate function values
    funValU(iteration) = .5 * norm(Xu - U_idx*Xu,'fro')^2 + lambda_1*norm(U_idx,1) +...
                         beta_1/2 * trace(U_idx'* L * U_idx); 

    if iteration > 1
        if funValU(iteration) < 1*10^-3
            break
        end
    end
    %fprintf('\nfunction value %.12f, iteration %d\n',funVal(iteration),iteration);
    if iteration > 100
        if funValU(iteration) < 1*10^-4 || funValU(iteration-1) == funValU(iteration) ...
                || funValU(iteration-1) - funValU(iteration) < 1*10^-4
            break;
        end
    end
    
end
end
