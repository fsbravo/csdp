% C-SDP relaxation for the quadratic assignment problem
%
% author: Jose Bravo-Ferreira
%         email: josesf@princeton.edu
%         Program in Applied and Computational Mathematics
%         Princeton University
%
%
% Note on usage: Function attempts to MINIMIZE trace(C*Q). If using for
% graph matching, use problem.A = A and problem.B = -B instead in order to
% maximize match instead.
%
% Input:
%   struct 'problem' with fields:
%     A      : first graph
%     B      : second graph (this should be the sparse graph)
%     P0     : starting permutation (default: identity)
%     maxits : maximum number of iterations
%     varsize: size of cliques used to create the PSD variables (default 2)
%     nprocs : number of processors to use in 'parallel' mode
%     stoch  : stochastic project (default 1)
%     mode   : admm (serial, default), parallel, ip (interior point, slow)
%
% Output:
%   struct 'results' with fields:
%     P                     : permutation matrix
%     P_stochastic          : permutation matrix with lowest cost
%     D                     : doubly stochastic matrix
%     optval                : objective value (lower bound)
%     upper_bound           : upper bound using P
%     upper_bound_stochastic: upper bound using P_stochastic
%     eta_vec               : vector of convergence criterion

function [results] = assign_csdp(problem)
    assert(isfield(problem, 'A'), 'Field "A" must be provided!');
    assert(isfield(problem, 'B'), 'Field "A" must be provided!');
    assert(size(problem.A,1)==size(problem.A,2), 'Graphs must be square!');
    assert(size(problem.B,1)==size(problem.B,2), 'Graphs must be square!');
    if ~isfield(problem, 'P0')
        fprintf('Initializing permutation to the identity.\n');
        problem.P0 = eye(size(problem.A,1));
    end
    if ~isfield(problem, 'maxits')
        warning('Field "maxits" not provided. Using maxits = 2000.');
        problem.maxits = 2000;
    end
    if ~isfield(problem, 'varsize')
        warning('Field "varsize" not provided. Using cliques of size 2.');
        problem.varsize = 2;
    end
    if ~isfield(problem, 'stoch')
        fprintf('Stochastic projection turned on.\n');
        problem.stoch = 1;
    end
    if ~isfield(problem, 'mode')
        warning('Field "mode" not provided. Using serial ADMM.');
        problem.mode = 'admm';
    else
        
    end
    % call appropriate solver
    if strcmp(problem.mode, 'admm')
        results = assign_csdp_admm(problem);
    elseif strcmp(problem.mode, 'parallel')
        if ~isfield(problem, 'nprocs')
            warning('Field "nprocs" not provided. Using serial admm solver instead.')
            results = assign_csdp_admm(problem);
        else
            results = assign_csdp_parallel(problem);
        end
    elseif strcmp(problem.mode, 'ip')
        warning('Interior point method is extremely slow for n>20! Use modes "admm" or "parallel" instead.')
        results = assign_csdp_ip(problem);
    else
        warning('Field "mode" not recognized! Using serial admm solver.');
        results = assign_csdp_admm(problem);
    end
end