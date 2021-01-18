function result = mosif_solve_with_clp(model, options, suppress_print)
%% Wrapper function to clp solver from OPTI distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%  model              - structure with formulation of optimization problem.
%                       For details, refer to function 'solve_optprob'.
%  options            - structure containing solver-specific options. Refer
%                       to the solver's manual for details.
%OUTPUTS:
%  result             - structure containing the following fields:
%   .x                - value of variable vector at solution
%   .f                - value of objective function at solution
%   .status           - character array showing exit condition. Values are:
%                        'optimal'
%                        'infeasible'
%                        'unbounded'
%                        'exceeded iterations'
%                        'other' (solver-specific reason)
%   .lm               - structure/vector with values of dual variables at 
%                       solution (currently not supported for all solvers)
%   .extra            - structure containing solver-specific information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% initialize output
result = struct('x', [], 'f', [], 'lm', [], 'status', 'unsuitable solver', 'extra', []);

%% ensure quadratic objective matrix is empty for LPs
Q = tril(model.Q);
if nnz(Q) == 0
    Q = [];
end

%% ensure linear constraints exist
A = model.A; 
bl = model.bl; 
bu = model.bu;
if isempty(A)
    A = sparse(1, length(model.x0));
    bl = -inf;
    bu = inf;
end

%% specify options if provided and suppress printing if need be
opt = options;
if suppress_print
    opt.display = 0;
end

%% solve problem
[result.x, result.f, exitflag, iter, result.lm] = clp(Q, model.q, ...
    A, bl, bu, model.xl, model.xu, opt);

%% record solution
result.extra = struct('status', exitflag, 'iter', iter);
if exitflag == 0
    result.status = 'optimal';
elseif exitflag == 1
    result.status = 'infeasible';
elseif exitflag == 2
    result.status = 'unbounded';
elseif exitflag == 3
    result.status = 'exceeded iterations';
else
    result.status = 'other';
end

end

