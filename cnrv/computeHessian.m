syms u11 u12 u21 u22 v11 v12 v13 v21 v22 v23;
U = [0 u12; u21 0];
%V = [v11 v12 v13;v21 v22 v23];
V = [0 0 0;0 0 0];
X = [0.128547373 0.887938318;0.166189151 0.888116686;0.212144432 0.897489234];

%VR = [(v12-0) (v13-v12); (0-v21) (v23-0)];

f = .5 * norm(X' - U*V,'fro')^2 + -1.00*norm(U,1); %+ 2.00*norm(V,1) + 3.00*l2l1norm(VR);
%u12,u21,v12,v13,v21,v23
H = hessian(f,[u12,u21])
%size(H)
eigen = eig(vpa(H));
ans = min(eigen)