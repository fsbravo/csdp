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

function [DD] = build_DD(n,BLKSZ)
    d = BLKSZ*n+1;
    DD = spalloc((BLKSZ*n)^2,d^2,(BLKSZ*n)^2);
    l = 1;
    for iii=1:BLKSZ*n
        for jjj=1:BLKSZ*n
            DD(l,(jjj-1)*d+iii) = 1;
            l = l+1;
        end
    end
end
