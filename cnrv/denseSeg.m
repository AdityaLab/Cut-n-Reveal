function t = denseSeg(s,width)
    % Condense the clusters into a single vector. 
    % @param s: a vector of dimension = num_timesteps X num_clusters. 
    % @param width: 
    
    [~, idx] = sort(sum(s,1),'descend');  %sort in order of the number of items in cluster. 
    s = s(:,idx);  %here idx stores the column indices of the original sum(s,1) matrix before sorting.
    
    for i=1:size(s,2)
        s(:,i) = i*s(:,i);
    end;
    
    s = sum(s,2);  %column wise sums.
    t = repmat(s,[1,width]);
end