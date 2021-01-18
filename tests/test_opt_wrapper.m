% test wrapper function for solvers
clear; clc;

%% user-defined options for testing
% names of solvers to test
solvers = {'ipopt', 'clp', 'gurobi'} ; %, 'mosek', 'snopt', 'minos', 'npsol', 'knitro'};

% optimization problem types to test
problems = {'model_LP_infeasible', 'model_LP_unbounded', 'model_LP_optimal', ...
    'model_QP_unconstrained', 'model_QP_infeasible', 'model_QP_optimal', ...
    'model_QCQP_infeasible', 'model_QCQP_optimal', 'model_NLP_infeasible', ...
    'model_NLP_optimal', 'model_MILP_infeasible', 'model_MILP_optimal', ...
    'model_MIQP_infeasible', 'model_MIQP_optimal', 'model_MIQCQP_infeasible', ...
    'model_MIQCQP_optimal'};

% printing to command window
suppress_print = 1;


%% prepare data
% get number of solvers and problem types
np = length(problems);
ns = length(solvers);

% get number of characters in each solver name (for nice printing later)
nchars = cellfun(@length, solvers);
nmax = max(nchars);

% get models and correct solutions for desired problems
models = cell(np, 1);
solutions = cell(np, 1);
for i = 1:np
    [models{i}, solutions{i}] = feval(problems{i});
end


%% test
% create table whose (i,j) entry shows how solver j passed test on problem i
Results = array2table(zeros(np, ns), 'VariableNames', solvers, 'RowNames', problems); 
% loop over all problem types
for i = 1:np
    fprintf('Testing solvers on problem type %s.\n', problems{i});
    
    % use all selected solvers to solve the given problem
    Results(i, :) = test_solvers_on_problem(models{i}, solvers, solutions{i}, suppress_print);
end


%% get aggregate information
fprintf('\nOverall statistics:\n');
for i = 1:ns
    solver_results = table2array(Results(:, i));
    solver_results(solver_results == -1) = [];
    fprintf('%s%s passed tests on %3.0f%% of suitable problems\n', ...
        upper(solvers{i}), blanks(nmax - nchars(i)), sum(solver_results) * 100 / length(solver_results));
end




%% test all given solvers on a problem instance
function pass = test_solvers_on_problem(model, solvers, sol, suppress_print)
ns = length(solvers);
pass = ones(1, ns);
opt.suppress_print = suppress_print;

% loop over solvers
for i = 1:ns
    % set solver
    opt.solver = solvers{i}; 
    
    % solve problem
    warning('off','all');
    solver_issues = false;
    try
        res = solve_optprob(model, opt);
    catch ME
        if strcmp(ME.message(1:6), 'Solver')
            solver_issues = true;
        end
    end
    warning('on','all');
    
    % check if problem type can be solved by the solver
    if solver_issues
        pass(i) = -1;
        continue;
    end
    
    % check solution status
    flag_status = strcmp(res.status, sol.status);
    % take care of special cases when solvers might not detect problem type
    if ~flag_status
        if strcmp(opt.solver, 'ipopt') && strcmp(sol.status, 'unbounded') ...
                && strcmp(res.status, 'exceeded iterations')
            flag_status = true;
        elseif strcmp(opt.solver, 'gurobi') && strcmp(res.extra.status, 'INF_OR_UNBD') ...
                && (strcmp(sol.status, 'infeasible') || strcmp(sol.status, 'unbounded'))
            flag_status = true;
        end
    end
    if ~flag_status
        fprintf('Solver %s failed solution status test\n', opt.solver);
    end
    
    % check variable vector and objective
    flag_variables = true;  
    flag_objective = true;
    if strcmp(sol.status, 'optimal')
        flag_variables = max(abs(sol.x - res.x)) < 1e-3;
        if ~flag_variables
            fprintf('Solver %s failed solution vector test\n', opt.solver);
        end
        flag_objective = abs(sol.f - res.f) < 1e-3;
        if ~flag_objective
            fprintf('Solver %s failed objective value test\n', opt.solver);
        end
        
    end
    
    % record information on passing the test
    pass(i) = flag_status & flag_variables & flag_objective;
end
pass = array2table(pass);
end





