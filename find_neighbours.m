function [neighbours] = find_neighbours(varlist,n)
    k = size(varlist,2);
    n_vars = size(varlist,1);
    neighbours = cell(1,k);
    idx = 1:n_vars;
    for i=1:k
        neighbours{i} = cell(1,n);
        for j=1:n
            neighbours{i}{j} = idx(varlist(:,i)==j);
        end
    end
end