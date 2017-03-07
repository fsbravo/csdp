%% matpart.m
% partitions matrix along dimension <dim> according to the the lengths
% specificied in <vec>. 
% E.g. given:
%   A = [1 2 3 4 5; 6 7 8 9 10]
%   dim = 2
%   vec = [2 2 1]
% produces:
%   part = cell(1,3)
% where: 
%   part{1} = [1 2; 6 7]
%   part{2} = [3 4; 8 9]
%   part{3} = [5; 10]
% NOTE: sum(vec) must be the same as size(A,dim)

function [part] = matpart(A,dim,vec)
    assert((dim==1)||(dim==2),'Dimension must be 1 or 2');
    assert(sum(vec)==size(A,dim),'Sum of block sizes must equal size of matrix');
    count = 1;
    for i=1:length(vec)
        val = vec(i);
        if dim==1
            part{i} = A(count:count+val-1,:);
        elseif dim==2
            part{i} = A(:,count:count+val-1);
        end
        count = count+val;
    end
end