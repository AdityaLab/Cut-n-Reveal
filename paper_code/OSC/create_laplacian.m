function[L] = create_laplacian(fileName, adjlist)
%py.importlib.import_module('networkx');
%counties= importdata(fileName,',');
%G = py.networkx.read_adjlist(adjlist, 'delimiter','\t');
%S = py.networkx.laplacian_matrix(G,counties.textdata(1,:));
%disp(S);
%[m,n] = S.shape;
%s = S.data;
%i = S.tocoo().row+1;
%j = S.indices+1;
L = load '../laplace.mat';
disp(L);
