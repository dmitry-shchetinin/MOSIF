function [result, extra] = solve_problem(model, options)
%% Wrapper function to various optimization solvers
% supported problem types
% NLP
%    minimize    f(x)
%    subject to  cl <= g(x) <= cu
%                bl <= A * x <= bu
%                xl <= x <= xu
% QP/LP
%    minimize    1/2 * x^T * Q * x + c^T * x
%    subject to  bl <= A * x <= bu
%                xl <= x <= xu
%
% supported and tested solvers: 
%      TOMLAB: snopt, knitro, minos, npsol
%      OPTI:   ipopt, clp
%      mosek
%
% INPUTS
%  model              - structure with formulation of optimization problem.
%                       Fields are (if not provided, it is set to empty or 
%                       default value):
%   .name             - problem name (character array)
%   .type             - problem type (character array). Supported types are
%                        'LP' - minimize (c^T * x + d) s.t. constraints
%                        'QP' - minimize (1/2 * x^T * Q * x + c^T * x) 
%                               s.t. constraints
%                        'NLP' (default) - minimize (f(x)) s.t. constraints
%   .Q                - quadratic coefficient matrix in objective (QPs)
%   .c                - linear coefficient vector in objective (LPs, QPs)
%   .funcs            - structure whose fields contain names of m-functions 
%                       that for NLP problems return the following values 
%                       at given point x:
%     .obj            - objective function (required for NLP)
%     .grad           - gradient of objective function
%     .con            - vector of nonlinear constraints
%     .jac            - Jacobian of nonlinear constraints
%     .h_obj          - Hessian of objective function
%     .h_con          - Hessian of nonlinear constraints
%   .x0               - starting point (REQUIRED)
%   .xl               - lower bounds on variables 
%   .xu               - upper bounds on variables
%   .A                - matrix of linear constraints
%   .bl               - lower bounds on linear constraints
%   .bu               - upper bounds on linear constraints
%   .cl               - lower bounds on nonlinear constraints
%   .cu               - upper bounds on nonlinear constraints
%   .Js               - sparsity pattern of nonlinear constraints' Jacobian
%                       (if not provided, it will be numerically estimated)
%   .Hs               - sparsity pattern of Hessian of Lagrangian (if not 
%                       provided, it will be numerically estimated)
%   .userdata         - structure containing data needed to compute values 
%                       of nonlinear functions
%  options            - structure containing various options for solving 
%                       the problem. The fields are:
%   .solver           - solver name. Supported and tested options include
%                        TOMLAB: 'snopt', 'knitro', 'minos', 'npsol'
%                        OPTI:   'ipopt' (default), 'clp'
%                        'mosek'
%   .printlevel       - level of command window printing ([] or 0 to 100,
%                       depending on the solver):
%                        [] - let each solver use its default (default)
%                        0 - nothing printed
%                        max - full info printed
%   .SOL              - row vector whose entries are parameter values of 
%                       SOL solvers (snopt, minos, npsol). Refer to solver 
%                       manuals for details
%   .KNITRO           - structure whose fields are options for TOMLAB's
%                       KNITRO. Refer to TOMLAB KNITRO manual for details
%   .IPOPT            - structure whose fields are options for IPOPT. Refer 
%                       to IPOPT manual for details
%   .CLP              - structure whose field are options for CLP. Refer to
%                       CLP manual for details
%   .MOSEK            - structure whose field are options for MOSEK. Refer 
%                       to MOSEK manual for details
% 
%OUTPUTS:
%  result             - structure containing the following fields:
%   .x                - value of variable vector at solution
%   .f                - value of objective function at solution
%   .status           - character array showing exit condition. Values are:
%                        'optimal'
%                        'infeasible'
%                        'unbounded'
%                        'exceeded iterations'
%                        'failed' (due to solver-specific reason)
%                        'incorrect solver'
%   .lm               - structure/vector with values of dual variables at 
%                       solution (currently not supported for all solvers)
%  extra              - structure containing solver-specific information
%
%DESCRIPTION OF INPUTS OF USER-PROVIDED NONLINEAR FUNCTIONS:
% User provided m-functions for computing objective, gradient, constraints, 
% Jacobian and Hessian of objective must have two inputs: X, STRUCT. 
% User provided m-function for computing Hessian of constraints must have 
% three inputs: X, LAMBDA, S. 
% Here, X is vector of variables, LAMBDA is vector of Lagrange multipliers 
% for nonlinear constraints, S is structure with field 'userdata', which 
% stores data needed to compute nonlinear functions. If no user data is 
% utilized by a function, S need not be touched.
% 
% Written by Dmitry Shchetinin (dmitry.v.shchetinin@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% preprocessing
% initialize outputs to empty
result = struct('x', [], 'f', [], 'status', 'failed', 'lm', []);
extra = [];

% check that required fields are present in model
check_inputs(model);
    
% fill out missing fields of model structure with empty or default values
model = populate_model(model);

% fill out missing fields of options structure with empty or default values
if nargin == 1
    options = struct();
end
options = populate_options(options);

% make sure all matrices are sparse
model.Q = sparse(model.Q);
model.A = sparse(model.A);
model.Js = sparse(model.Js);
model.Hs = sparse(model.Hs);

% check if problem belongs to class of standard convex problems
if strcmp(model.type,'QP') || strcmp(model.type,'LP')
    isCP = true;
else
    isCP = false;
end
model.isCP = isCP;

% check if chosen solver can solve a given problem
if ~isCP && ~strcmp(options.solver,'ipopt') && ~strcmp(options.solver,'snopt') && ...
        ~strcmp(options.solver,'minos') && ~strcmp(options.solver,'npsol') && ...
        ~strcmp(options.solver,'knitro')
    warning('NLP problem cannot be solved by %s. Select a different solver.',options.solver);
    result.status = 'incorrect solver';
    return;
end

% make sure that all coefficient vectors are column vectors
model.c = make_col_vec(model.c);
model.bl = make_col_vec(model.bl);
model.bu = make_col_vec(model.bu);
model.cl = make_col_vec(model.cl);
model.cu = make_col_vec(model.cu);
model.xl = make_col_vec(model.xl);
model.xu = make_col_vec(model.xu);

% extract model dimensions
nv = numel(model.x0); % number of variables
nl = size(model.A,1); % number of linear constraints
nn = 0;               % number of nonlinear constraints
if strcmp(model.type,'NLP') && ~isempty(model.funcs.con)
    if isempty(model.cl) && isempty(model.cu)
        nn = numel(feval(model.funcs.con, model.x0, model));
    else
        nn = max(numel(model.cl), numel(model.cu));
    end
end

% estimate Jacobian sparsity pattern if needed
Js = model.Js;
if strcmp(model.type,'NLP') && isempty(Js) && ~isempty(model.funcs.con)
    Js = Jacobian_sparsity_pattern(model);
end

% estimate Hessian sparsity pattern if needed
if isCP
    Hs = model.Q;
    Hs(Hs ~= 0) = 1;
elseif ~isempty(model.Hs)
    Hs = model.Hs;
elseif strcmp(model.type,'NLP') && ~isempty(model.funcs.h_con)
    Hs = Hessian_sparsity_pattern(model);
else 
    Hs = [];
end

% set infinite bounds on linear constraints if not provided
model = populate_field(model,'bl', -inf(nl,1));
model = populate_field(model,'bu', inf(nl,1));

% set infinite bounds on nonlinear constraints if not provided
model = populate_field(model,'cl', -inf(nn,1));
model = populate_field(model,'cu', inf(nn,1));



%% solve problem
if strcmp(options.solver,'ipopt')
    %% solve problem by IPOPT (default solver)
    % set infinite bounds on variables if not provided
    model = populate_field(model,'xl', -inf(nv,1));
    model = populate_field(model,'xu', inf(nv,1));
  
    % create model structure required by ipopt
    auxdata = struct('userdata', model.userdata, 'funcs', model.funcs, ...
        'A', model.A, 'Q', model.Q, 'c', model.c, 'isCP', model.isCP);
    opt = struct('lb', model.xl, 'ub', model.xu, 'cl', [model.cl; model.bl], ...
        'cu', [model.cu; model.bu]);
    
    % specify options if provided
    if ~isempty(options.IPOPT)
        opt.ipopt = options.IPOPT;
    end
    if ~isempty(options.printlevel) && (isempty(options.IPOPT) || ...
            ~isfield(opt.ipopt, 'print_level'))
        opt.ipopt.print_level = options.printlevel;
    end
    
    % specify function handles for objective function and its gradient
    funcs.objective = @(x) objective(x, auxdata);
    funcs.gradient  = @(x) gradient(x, auxdata);
    
    % specify function handles for constraints if exist, and their Jacobian
    if ~isempty(model.funcs.con) || ~isempty(model.A)
        funcs.constraints = @(x) constraints(x, auxdata);
        funcs.jacobian    = @(x) jacobian(x, auxdata);
        
        % provide Jacobian sparsity pattern
        Js = [Js; model.A];
        Js(Js ~= 0) = 1;
        funcs.jacobianstructure = @(d) (Js);
    end
    
    % specify function handles for Hessian if provided
    if ~isempty(Hs) && (~isfield(options.IPOPT, 'hessian_approximation') || ...
            ~strcmp(options.IPOPT.hessian_approximation, 'limited-memory'))
        funcs.hessian  = @(x, sigma, lambda) hessian(x, sigma, lambda, auxdata);
        
        % provide Hessian sparsity pattern
        funcs.hessianstructure = @(d) (tril(Hs));
    else
        opt.ipopt.hessian_approximation = 'limited-memory';
    end
    
    % solve the problem
    [result.x, extra] = ipopt(model.x0, funcs, opt);
    exitflag = extra.status;
    
    % record results
    result.f = funcs.objective(result.x);
    result.lm = struct('zl', extra.zl, 'zu', extra.zu, 'lambda', extra.lambda);
    result.extra = extra;
    if exitflag == 0 || exitflag == 1
        result.status = 'optimal';
    elseif exitflag == -1
        result.status = 'exceeded iterations';
    elseif exitflag == 2
        result.status = 'infeasible';
    elseif exitflag == 3 || exitflag == 4 || exitflag == 5 || ...
            exitflag == -2 || exitflag == -3 || exitflag == -10
        result.status = 'unbounded';
    end
    
elseif strcmp(options.solver, 'clp') && isCP
    %% solve problem by CLP
    % esnure quadratic objective matrix is empty for LPs
    Q = tril(model.Q);
    if numel(Q) == 0
        Q = [];
    end
    
    % ensure linear constraints exist
    A = model.A; bl = model.bl; bu = model.bu;
    if isempty(A)
        A = sparse(1,nv);
        bl = -inf;
        bu = inf;
    end
    
    % specify options if provided
    opt = options.CLP;
    if ~isempty(options.printlevel) && (isempty(opt) || ~isfield(opt, 'display'))
        opt.display = options.printlevel;
    end
    
    % solve problem
    [result.x, result.f, exitflag, iter, result.lm] = clp(Q, model.c, ...
        A, bl, bu, model.xl, model.xu, opt);
    
    % record solution
    extra = struct('status', exitflag, 'iter', iter);
    if exitflag == 0
        result.status = 'optimal';
    elseif exitflag == 1
        result.status = 'infeasible';
    elseif exitflag == 2
        result.status = 'unbounded';
    elseif exitflag == 3
        result.status = 'exceeded iterations';
    end
    
elseif strcmp(options.solver, 'mosek') && isCP
    %% solve problem by MOSEK
    % define problem structure in mosek format
    Prob = struct('c', model.c, 'blc', model.bl, 'buc', model.bu, ...
        'blx', model.xl, 'bux', model.xu);
    
    % add quadratic objective matrix
    if numel(model.Q) > 0
        [Prob.qosubi,Prob.qosubj,Prob.qoval] = find(tril(model.Q));
    end

    % add linear constraints matrix
    Prob.a = model.A;
    if isempty(Prob.a)
        Prob.a = sparse(0,nv);
    end
    
    % add starting point
    Prob.sol.itr.xx = model.x0;
    
    % set print level to command window
    if isempty(options.printlevel)
        cmd = sprintf('minimize info statuskeys(1) symbcon');
    else
        cmd = sprintf('minimize info echo(%d) statuskeys(1) symbcon',options.printlevel);
    end
    
    % solve problem
    [rcode, extra] = mosekopt(cmd, Prob, options.MOSEK);
    
    % extract solution
    hasSol = isfield(extra,'sol');
    if hasSol
        result.x = extra.sol.itr.xx;
        result.f = extra.sol.itr.pobjval;
    end
    
    % get symbolic constants
    if isfield(extra,'symbcon')
        sc = extra.symbcon;
    else
        [~,r] = mosekopt('symbcon');
        sc = r.symbcon;
    end

    % extract exit information
    if hasSol
        solsta = extra.sol.itr.solsta;
        prosta = extra.sol.itr.prosta;
        if rcode == 0 && solsta == sc.MSK_SOL_STA_OPTIMAL
            result.status = 'optimal';
        elseif prosta == sc.MSK_PRO_STA_PRIM_INFEAS
            result.status = 'infeasible';
        elseif prosta == sc.MSK_PRO_STA_DUAL_INFEAS
            result.status = 'unbounded';
        end
    elseif rcode == sc.MSK_RES_TRM_MAX_ITERATIONS
        result.status = 'exceeded iterations';
    end
    
else
    %% solve problem by Tomlab
    % create problem structure
    sigma = 1;
    if isCP
        Prob = conAssign(@(x) objective(x, model), @(x) gradient(x, model), ...
            @(x, lambda) hessian(x, sigma, lambda, model), [], model.xl, ...
            model.xu, model.name, model.x0, [], [], model.A, model.bl, ...
            model.bu, [], [], [], Js, model.cl, model.cu);
    else
        Prob = conAssign(model.funcs.obj, model.funcs.grad, model.funcs.h_obj, ...
            [], model.xl, model.xu, model.name, model.x0, [], [], model.A, ...
            model.bl, model.bu, model.funcs.con, model.funcs.jac, ...
            model.funcs.h_con, Js, model.cl, model.cu);
    end
    Prob.LargeScale = 1; % by default assume large sparse problem
    Prob.d2LPattern = Hs;
    
    % provide extra data to pass to functions
    Prob.userdata=model.userdata;
    
    % specify solver options if provided
    if ~isempty(options.SOL) && (strcmp(options.solver,'snopt') || ...
        strcmp(options.solver,'minos') || strcmp(options.solver,'npsol'))
        Prob.SOL.optPar = options.SOL;
    elseif ~isempty(options.KNITRO) && strcmp(options.solver,'knitro') 
        Prob.KNITRO.options = options.KNITRO;
    end
    
    % solve the problem
    extra = tomRun(options.solver, Prob, options.printlevel);
    exitflag = extra.ExitFlag;
    
    % record result
    result.x = extra.x_k;
    result.f = extra.f_k;
    result.lm = extra.v_k;
    if exitflag == 4
        result.status = 'infeasible';
    elseif exitflag == 2 || result.f < -1e10
        result.status = 'unbounded';
    elseif exitflag == 1
        result.status = 'exceeded iterations';
    elseif exitflag == 0 || (strcmp(options.solver,'npsol') && exitflag == 3)
        result.status = 'optimal';
    end
end

end



%% preprocessing functions
% add model fields
function model = populate_model(model)
    n = numel(model.x0);
    model = populate_field(model,'Q',sparse(n,n));
    model = populate_field(model,'c',zeros(n,1));
    model = populate_field(model,'A',[]);
    model = populate_field(model,'bl',[]);
    model = populate_field(model,'bu',[]);
    model = populate_field(model,'xl',[]);
    model = populate_field(model,'xu',[]);
    model = populate_field(model,'x0',[]);
    model = populate_field(model,'cl',[]);
    model = populate_field(model,'cu',[]);
    model = populate_field(model,'Js',[]);
    model = populate_field(model,'Hs',[]);
    model = populate_field(model,'funcs',[]);
    model = populate_field(model,'name','problem');
    model = populate_field(model,'type','NLP');
    model = populate_field(model,'userdata',[]);
    model.funcs = populate_field(model.funcs,'obj',[]);
    model.funcs = populate_field(model.funcs,'grad',[]);
    model.funcs = populate_field(model.funcs,'h_obj',[]);
    model.funcs = populate_field(model.funcs,'con',[]);
    model.funcs = populate_field(model.funcs,'jac',[]);
    model.funcs = populate_field(model.funcs,'h_con',[]);
end

% add option fields
function options = populate_options(options)
    options = populate_field(options,'solver','ipopt');
    options = populate_field(options,'printlevel',[]);
    options = populate_field(options,'SOL',[]);
    options = populate_field(options,'KNITRO',[]);
    options = populate_field(options,'IPOPT',[]);
    options = populate_field(options,'CLP',[]);
    options = populate_field(options,'MOSEK',[]);
end

% check that required fields are present in model
function check_inputs(model)
if ~isfield(model,'x0') || isempty(model.x0)
    error('Starting point must be provided.');
end
end

% populate given field of a structure
function s = populate_field(s,fieldname, fieldvalue)
    if ~isfield(s,fieldname) || isempty(s.(fieldname))
        s.(fieldname)=fieldvalue;
    end
end

% ensure that vector is a column vector
function v = make_col_vec(v)
    if ~iscolumn(v)
        v = transpose(v);
    end
end


%% sparsity patterns
function Js = Jacobian_sparsity_pattern(model)
    % get number of variables
    n = numel(model.x0);
    % compute sum of absolute values of Jacobian for 3 random values of x
    Js = abs(feval(model.funcs.jac, rand(n,1), model)) + ...
        abs(feval(model.funcs.jac, rand(n,1), model)) + ...
        abs(feval(model.funcs.jac, rand(n,1), model));
    Js(Js ~= 0) = 1;
    Js = sparse(Js);
end

function Hs = Hessian_sparsity_pattern(model)
    % get number of variables
    nv = numel(model.x0);
    % compute sum of absolute values of Hessian for 3 random values of x
    Hs = abs(feval(model.funcs.h_obj, rand(nv,1), model))+...
        abs(feval(model.funcs.h_obj, rand(nv,1), model))+...
        abs(feval(model.funcs.h_obj, rand(nv,1), model));
    if ~isempty(model.funcs.h_con)
        % get number of nonlinear constraints
        nn = numel(feval(model.funcs.con, model.x0, model));
        % compute sum of absolute values of Hessian for 3 random values of x
        Hs=Hs+abs(feval(model.funcs.h_con, rand(nv,1), rand(nn,1), model)) + ...
            abs(feval(model.funcs.h_con, rand(nv,1), rand(nn,1), model)) + ...
            abs(feval(model.funcs.h_con, rand(nv,1), rand(nn,1), model));
    end
    Hs(Hs ~= 0) = 1;
    Hs = sparse(Hs);
end


%% callback functions for computing objective, constraints, and derivatives
function f = objective(x, auxdata)
    if auxdata.isCP
        f = (x' * auxdata.Q * x) / 2 + auxdata.c' * x;
    else
        f = feval(auxdata.funcs.obj, x, auxdata);
    end
end

function df = gradient(x, auxdata)
    if auxdata.isCP
        df = (auxdata.Q + auxdata.Q') * x / 2 + auxdata.c;
    else
        df = feval(auxdata.funcs.grad, x, auxdata);
    end
end

function c = constraints(x, auxdata)
    if ~auxdata.isCP && ~isempty(auxdata.funcs.con)
        cnonlin = feval(auxdata.funcs.con, x, auxdata);
    else
        cnonlin = [];
    end
    if ~isempty(auxdata.A)
        clin = auxdata.A * x;
    else
        clin = [];
    end
    c = [cnonlin; clin];
end

function J = jacobian(x, auxdata)
    if ~auxdata.isCP && ~isempty(auxdata.funcs.jac)
        J = [feval(auxdata.funcs.jac, x, auxdata); auxdata.A];
    else
        J = auxdata.A;
    end
end

function H = hessian(x, sigma, lambda, auxdata)
    if auxdata.isCP
        H = sigma * (auxdata.Q + auxdata.Q') / 2;
    else
        if ~isempty(auxdata.A)
            nl = size(auxdata.A,1); % number of linear constraints
            lambda = lambda(1:end - nl);
        end
        H = sigma * feval(auxdata.funcs.h_obj, x, auxdata);
        if ~isempty(auxdata.funcs.h_con)
            H = H + feval(auxdata.funcs.h_con, x, lambda, auxdata);
        end
    end
    H = tril(H);
end