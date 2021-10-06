function L = l2l1norm(x)
% L2L1 norm

    L = 0;
    for i=1:size(x,2)  %iterate over each column, calculate the 2-norm and add to L hence calculating the L1 norm.
        L = L + norm(x(:,i));
    end
end
