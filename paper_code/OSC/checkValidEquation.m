function [fV]= checkValidEquation(A,B,C)
fV=0;
vA = nnz(isnan(A))+nnz(isinf(A))+all(A(:)==0);
vB = nnz(isnan(B))+nnz(isinf(B))+all(B(:)==0);
vC = nnz(isnan(C))+nnz(isinf(C))+all(C(:)==0);
if ~vA && ~vB && ~vC
	eA = eig(A);
	eB= eig(B);
	%common = nnz(ismember(abs(eA),abs(eB)));
	ttl=0;
	for i=1:size(eA,2)
		for j=1:size(eB,2)
			if (eA(i) + eB(j)) < 1*10^-10
				ttl=ttl+1;
				break;
			end
		end
		if ttl>0
			break;
		end
	end
	if ttl==0
		fV=1;
	end
end
end