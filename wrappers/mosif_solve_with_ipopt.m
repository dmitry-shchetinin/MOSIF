function result = mosif_solve_with_ipopt(model, options, suppress_print)
%% Wrapper function to ipopt solver from OPTI distribution
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

%% create model structure required by ipopt
auxdata = struct('userdata', model.userdata, 'funcs', model.funcs, ...
    'A', model.A, 'Q', model.Q, 'q', model.q, 'p', model.p, 'isNLP', ...
    model.isNLP, 'quadcon', model.quadcon);
opt = struct('lb', model.xl, 'ub', model.xu, 'cl', [model.cl; model.bl], ...
    'cu', [model.cu; model.bu]);

%% specify options if provided
if ~isempty(options)
    opt.ipopt = options;
end

%% suppress output to command window if desired
if suppress_print
    opt.ipopt.print_level = 0;
end

%% specify function handles for objective function and its gradient
funcs.objective = @(x) mosif_objective(x, auxdata);
funcs.gradient  = @(x) mosif_gradient(x, auxdata);

%% specify function handles for constraints if provided, and their Jacobian
if ~isempty(model.funcs.con) || ~isempty(model.A)
    funcs.constraints = @(x) mosif_constraints(x, auxdata);
    funcs.jacobian    = @(x) mosif_jacobian(x, auxdata);

    % provide Jacobian sparsity pattern
    Js = [model.Js; mosif_set_entries_to_ones(model.A)];
    funcs.jacobianstructure = @(d) (Js);
end

%% specify function handles for Hessian if provided
if ~isempty(model.Hs) && (isempty(options) || ~isfield(options, 'hessian_approximation') || ...
        ~strcmp(options.hessian_approximation, 'limited-memory'))
    funcs.hessian  = @(x, sigma, lambda) mosif_hessian(x, sigma, lambda, auxdata);

    % provide Hessian sparsity pattern
    funcs.hessianstructure = @(d) (tril(model.Hs));
else
    opt.ipopt.hessian_approximation = 'limited-memory';
end

%% solve the problem
[result.x, extra] = ipopt(model.x0, funcs, opt);

%% record results
result.f = funcs.objective(result.x);
result.lm = struct('zl', extra.zl, 'zu', extra.zu, 'lambda', extra.lambda);
result.extra = extra;
% record status of solution
exitflag = extra.status;
if exitflag == 0 || exitflag == 1
    result.status = 'optimal';
elseif exitflag == -1
    result.status = 'exceeded iterations';
elseif exitflag == 2
    result.status = 'infeasible';
elseif exitflag == 3 || exitflag == 4 || exitflag == 5 || ...
        exitflag == -2 || exitflag == -3 || exitflag == -10
    result.status = 'unbounded';
else
    result.status = 'other';
end

end

