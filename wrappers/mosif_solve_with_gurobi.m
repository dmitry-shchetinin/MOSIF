function result = mosif_solve_with_gurobi(model, options, suppress_print)
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

%% convert conic constraints to quadratic
model = mosif_cone2quadratic(model);

%% initialize gurobi model structure
gurmod = struct();

%% populate fields describing objective function
gurmod.Q = model.Q;
gurmod.obj = model.q;
gurmod.objcon = model.p;

%% populate fields describing variable bounds and types
gurmod.lb = model.xl;
gurmod.ub = model.xu;
gurmod.vtype = model.xtype;

%% populate fields describing linear constraints
% equality constraints
if ~isempty(model.A)
    idx = model.bl == model.bu;
    A = model.A(idx, :);
    rhs = model.bl(idx);
    sense = repelem('=', sum(idx));
    % inequality constraints with upper bounds
    idx = model.bu ~= inf;
    A = [A; model.A(idx, :)];
    rhs = [rhs; model.bu(idx)];
    sense = [sense, repelem('<', sum(idx))];
    % inequality constraints with upper bounds
    idx = model.bl ~= -inf;
    A = [A; model.A(idx, :)];
    rhs = [rhs; model.bl(idx)];
    sense = [sense, repelem('>', sum(idx))];
    % record into model
    gurmod.A = A;
    gurmod.rhs = rhs;
    gurmod.sense = sense;
else
    gurmod.A = sparse(1,length(model.x0));
end

%% populate fields describing quadratic constraints
n_quadcon = length(model.quadcon);
for i = 1:n_quadcon
    gurmod.quadcon(i).Qc = model.quadcon(i).Qc;
    gurmod.quadcon(i).q = model.quadcon(i).qc;
    gurmod.quadcon(i).rhs = model.quadcon(i).rhs;
end

%% suppress printing if need be
if suppress_print
    options.OutputFlag = 0;
end

%% solve the problem
gurres = gurobi(gurmod, options);

%% record results
if isfield(gurres, 'objval')
    result.f = gurres.objval;
end
if isfield(gurres, 'x')
    result.x = gurres.x;
end
result.extra = gurres;
% record status of solution
if strcmp(gurres.status, 'OPTIMAL') || strcmp(gurres.status, 'SUBOPTIMAL')
    result.status = 'optimal';
elseif strcmp(gurres.status, 'ITERATION_LIMIT')
    result.status = 'exceeded iterations';
elseif strcmp(gurres.status, 'INFEASIBLE')
    result.status = 'infeasible';
elseif strcmp(gurres.status, 'UNBOUNDED')
    result.status = 'unbounded';
else
    result.status = 'other';
end

end

