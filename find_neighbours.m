% C-SDP: SDP relaxation techniques for the Quadratic Assignment Problem
% Copyright (C) 2017  Jose Bravo Ferreira
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
