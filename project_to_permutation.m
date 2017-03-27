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

% add white noise to a doubly stochastic matrix before using hungarian
% algorithm to obtain better permutation matrix
function [P] = project_to_permutation(D,A,B)
    fprintf('Projecting to a permutation\n');
    n = size(A,1);
    itsperstep = 1000;
    
    max_noise = 0.5;
    n_noise = 10;
    noise_v = linspace(0,max_noise,n_noise+1);
    noise_v = noise_v(2:end);
    
    P = closest_permutation(D);
    best_val = trace(A*P*B*P');
    for i=1:n_noise
        noise = noise_v(i);
        for j=1:itsperstep
            D_n = D + normrnd(0,noise,n,n);
            P_n = closest_permutation(D_n);
            obj = trace(A*P_n*B*P_n');
            if obj<best_val
                best_val = obj;
                P = P_n;
            end
        end
    end
end
