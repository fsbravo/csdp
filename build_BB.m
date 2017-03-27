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

function [BB_col] = build_BB(n,BLKSZ)
    BB_col = cell(1,BLKSZ);
    d = BLKSZ*n+1;
    for j=1:BLKSZ
        BB = spalloc(n,d^2,2*n);
        for iii=1:n
            temp = zeros(d); temp(iii+(j-1)*n,d) = 1/2; temp(d,iii+(j-1)*n) = 1/2;
            BB(iii,:) = temp(:)';
        end
        BB_col{j} = BB;
    end
end
