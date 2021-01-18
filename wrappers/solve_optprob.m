function result = solve_optprob(model, options)
%% Wrapper function to various optimization solvers
% supported problem types
% NLP
%    minimize    f(x)
%    subject to  cl <= g(x) <= cu
%                bl <= A * x <= bu
%                xl <= x <= xu
% QCQP/QP/LP + SOCP
%    minimize    x^T * Q * x + q^T * x + p
%    subject to  bl <= A * x <= bu
%                x^T * Qc_i * x + qc_i^T * x <= rhs_i (for all i)
%                x_i^2 + ... x_j^2 <= x_k^2  (cones)
%                x_i^2 + ... + x_j^2 <= x_k * x_m (rotated cones)
%                xl <= x <= xu
%
% supported and tested solvers: 
%      TOMLAB: snopt, knitro, minos, npsol
%      OPTI:   ipopt, clp
%      mosek
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%  model              - structure with formulation of optimization problem.
%                       Fields are (if not provided, it is set to empty or 
%                       default value):
%   .name             - problem name (character array)
%   .Q                - sparse quadratic coefficient matrix in objective
%   .q                - linear coefficient vector in objective
%   .p                - constant term in the objective function
%   .funcs            - structure whose fields contain names of m-functions 
%                       that for NLP problems return the following values 
%                       at given point x:
%     .obj            - objective function (required for NLP)
%     .grad           - gradient of objective function
%     .con            - vector of nonlinear constraints
%     .jac            - Jacobian of nonlinear constraints
%     .hess_obj       - Hessian of objective function
%     .hess_con       - Hessian of nonlinear constraints
%   .x0               - starting point (REQUIRED)
%   .xl               - lower bounds on variables 
%   .xu               - upper bounds on variables
%   .xtype            - character array of variable types. Elements can be:
%                       'C' - continuous (default)
%                       'B' - binary
%                       'I' - integer
%   .A                - sparse matrix of linear constraints
%   .bl               - lower bounds on linear constraints
%   .bu               - upper bounds on linear constraints
%   .cl               - lower bounds on nonlinear constraints
%   .cu               - upper bounds on nonlinear constraints
%   .Js               - sparsity pattern of nonlinear constraints' Jacobian
%                       (if not provided, it will be numerically estimated)
%   .Hs               - sparsity pattern of Hessian of Lagrangian (if not 
%                       provided, it will be numerically estimated)
%   .quadcon          - structure array representing quadratic constraints.
%                       Each element models a single constraint and is a 
%                       structure with the following fields:
%     .Qc             - sparse quadratic coefficient matrix
%     .qc             - sparse linear coefficient vector
%     .rhs            - right-hand side scalar
%   .cones            - structure array representing second-order cone 
%                       constraints. Each entry models a single constraint 
%                       and is a structure with the following fields:
%     .idx_le         - vector of indices of variables that belong to the 
%                       "less or equal" side of the constraint
%     .idx_ge         - vector of indices of variables that belong to the 
%                       "greater or equal" side of the constraint. It must
%                       have a single element for a cone constraint and two
%                       elements for a rotated cone constraint.
%   .userdata         - structure containing data needed to compute values 
%                       of nonlinear functions
%  opts               - structure containing various general and solver-
%                       specific options for solving the optimization 
%                       problem. If it is not provided or is empty, all 
%                       fields are initialized to their defaults. If a 
%                       field is not provided by the user, it will be 
%                       initizalized to its default. The fields are 
%                       (default values are given after field names):
%   .solver   'ipopt' - solver name. Supported and tested options are:
%                       TOMLAB: 'snopt', 'knitro', 'minos', 'npsol'
%                       OPTI:   'ipopt' , 'clp'
%                       'mosek'
%   .suppress_print 0 - decide on what to print to command window:
%                       0 - solver will use its default printing level
%                       1 - nothing is printed (NOTE: if selected, this
%                           will overwrite printing level that the user
%                           might have set in the solver options).
%   .ipopt        []  - structure whose fields are options for IPOPT. 
%   .clp          []  - structure whose fields are options for CLP. 
%   .mosek        []  - structure whose fields are options for MOSEK. 
%   .tomlab       []  - structure whose fields are options for solvers 
%                       that are part of Tomlab distribution. Fields are:
%     .sol        []  - row vector whose entries are options of Tomlab's 
%                       SOL solvers (snopt, minos, npsol). 
%     .knitro     []  - structure whose fields are options for Tomlab's 
%                       KNITRO.
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


%% fill out missing fields of options structure with empty/default values
if nargin == 1
    options = struct();
end
options = mosif_populate_options(options);

%% preprocess model to make sure all fields are consistent and correct
model = mosif_preprocess_model(model, options.solver);

%% solve problem
if strcmp(options.solver, 'ipopt')
    result = mosif_solve_with_ipopt(model, options.ipopt, options.suppress_print);
elseif strcmp(options.solver, 'clp')
    result = mosif_solve_with_clp(model, options.clp, options.suppress_print);
elseif strcmp(options.solver, 'gurobi')
    result = mosif_solve_with_gurobi(model, options.gurobi, options.suppress_print);    
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