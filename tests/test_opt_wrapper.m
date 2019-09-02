% test wrapper function for solvers
clear; clc;

%% user-defined options for testing
% names of solvers to test
solvers = {'ipopt', 'clp', 'mosek', 'snopt', 'minos', 'npsol', 'knitro'};

% optimization problem types to test
problems = {'LP_infeasible', 'LP_unbounded', 'LP_optimal', 'QP_unconstrained', ...
    'QP_infeasible', 'QP_optimal', 'NLP_infeasible', 'NLP_optimal'};

% printing to command window
printlevel = 0;


%% prepare data
% get number of solvers and problem types
np = length(problems);
ns = length(solvers);

% get number of characters in each solver name (for nice printing later)
nchars = cellfun(@length,solvers);
nmax = max(nchars);

% get testcases for desired problem types
models = cell(np, 1);
solutions = cell(np, 1);
for i = 1:np
    [models{i}, solutions{i}] = feval(problems{i});
end


%% test
Results = zeros(ns, np); % matrix whose (i,j) entry shows how solver i passed test on problem j
% loop over all problem types
for i = 1:np
    fprintf('Testing solvers on problem type %s.\n', problems{i});
    
    if i ==6
        rrr=5;
    end
    
    % check all solvers
    Results(:,i) = test_solvers_on_problem(models{i}, solvers, solutions{i}, printlevel);
end

% compute aggregate information
fprintf('\nOverall statistics:\n');
for i = 1:ns
    fprintf('%s%s passed tests on %3.0f%% of problems\n', ...
        upper(solvers{i}), blanks(nmax - nchars(i)), sum(Results(i,:)) * 100 / np);
end




%% test all given solvers on a problem instance
function pass = test_solvers_on_problem(model, solvers, sol, printlevel)
ns = length(solvers);
pass = true(ns, 1);
opt.printlevel = printlevel;

% loop over solvers
for i = 1:ns
    % set solver
    opt.solver = solvers{i}; 
    
    % solve problem
    warning('off','all');
    res = solve_problem(model, opt);
    warning('on','all');
    
    % check if problem type can be solved by the solver
    if strcmp(res.status, 'incorrect solver')
        continue;
    end
    
    % check solution status
    flag1 = strcmp(res.status, sol.status);
    if ~flag1
        fprintf('Solver %s failed solution status test\n', opt.solver);
    end
    
    % check variable vector and objective
    flag2 = true;  flag3 = true;
    if strcmp(sol.status, 'optimal')
        flag2 = max(abs(sol.x - res.x)) < 1e-5;
        if ~flag2
            fprintf('Solver %s failed solution vector test\n', opt.solver);
        end
        flag3 = abs(sol.f - res.f) < 1e-5;
        if ~flag3
            fprintf('Solver %s failed objective value test\n', opt.solver);
        end
        
    end
    
    % record information on passing the test
    pass(i) = flag1 & flag2 & flag3;
end
end



%% LP problems
function [model, sol] = LP_infeasible()
% define problem
model.type = 'LP';
model.name = 'LP_infeas';
model.x0 = zeros(4,1);
model.c = ones(4,1);
model.xu = zeros(4,1);
model.A = sparse(ones(1,4));
model.bl = 1;
model.bu = 1;
% define correct solution
sol.status = 'infeasible';
sol.x = [];
sol.f = [];
end

function [model, sol] = LP_unbounded()
% define problem
model.type = 'LP';
model.name = 'LP_unbnd';
model.x0 = zeros(4,1);
model.c = ones(4,1);
model.xu = ones(4,1);
model.A = sparse(ones(1,4));
model.bu = 0;
% define correct solution
sol.status = 'unbounded';
sol.x = [];
sol.f = [];
end

function [model, sol] = LP_optimal()
% define problem
model.type = 'LP';
model.name = 'LP_optml';
model.x0 = zeros(4,1);
model.c = [1; 2; 3; 4];
model.xl = zeros(4,1);
model.A = sparse(ones(1,4));
model.bl = 1;
model.bu = 1;
% define correct solution
sol.status = 'optimal';
sol.x = [1; 0; 0; 0];
sol.f = 1;
end


%% QP problems
function [model, sol] = QP_unconstrained()
% define problem
model.type = 'QP';
model.name = 'QP_uncnstd';
model.x0 = [1; 2; 3; 4];
model.Q = speye(4);
% define correct solution
sol.status = 'optimal';
sol.x = zeros(4,1);
sol.f = 0;
end

function [model, sol] = QP_infeasible()
% define problem
model.type = 'QP';
model.name = 'QP_infeas';
model.x0 = zeros(4,1);
model.Q = speye(4);
model.c = ones(4,1);
model.xu = zeros(4,1);
model.A = sparse(ones(1,4));
model.bl = 1;
model.bu = 1;
% define correct solution
sol.status = 'infeasible';
sol.x = [];
sol.f = [];
end

% function [model, sol] = QP_unbounded()   % NOTE: clp fails on this
% % define problem
% model.type = 'QP';
% model.name = 'QP_unbnd';
% model.x0 = zeros(4,1);
% model.Q = sparse(1,1,1,4,4);
% model.c = [0; 1; 1; 1];
% model.A = sparse([1, 0, 0, 0]);
% model.bl = 1;
% model.bu = 1;
% % define correct solution
% sol.status = 'unbounded';
% sol.x = [];
% sol.f = [];
% end

function [model, sol] = QP_optimal()
% define problem
model.type = 'QP';
model.name = 'QP_optml';
model.x0 = zeros(4,1);
model.Q = sparse(1,1,1,4,4);
model.c = [0; 1; 1; 1];
model.A = sparse([0, 0, 1, -2]);
model.bl = 0;
model.bu = 0;
model.xl = -ones(4,1);
model.xu = ones(4,1);
% define correct solution
sol.status = 'optimal';
sol.x = [0; -1; -1; -0.5];
sol.f = -2.5;
end


%% NLP problems
%the following problem is solved (bounds are changed to create dif. problem types):
% min exp(x1) * (4*x1^2 + 2*x2^2 + 4*x1*x2 + 2*x2 + 1)
% s.t. -10<=x1<=10
%      -10<=x2<=10
%      x1+x2=0
%      -x1*x2 + x1 + x2>=1.5
%      x1*x2>=-10
function [model, sol] = NLP_infeasible()
% define problem
model.type = 'NLP';
model.name = 'NLP_infeas';
model.x0 = zeros(2,1);
model.funcs.obj = @prob_obj;
model.funcs.grad = @prob_grad;
model.funcs.h_obj = @prob_hobj;
model.funcs.con = @prob_con;
model.funcs.jac = @prob_jac;
model.funcs.h_con = @prob_hcon;
model.xu = [-10; -10];
model.A = sparse([1, 1]);
model.bl = 0;
model.bu = 0;
model.cl = [1.5; -10];
model.cu = [inf; inf];
% define correct solution
sol.status = 'infeasible';
sol.x = [];
sol.f = [];
end

function [model, sol] = NLP_optimal()
% define problem
model.type = 'NLP';
model.name = 'NLP_optml';
model.x0 = zeros(2,1);
model.funcs.obj = @prob_obj;
model.funcs.grad = @prob_grad;
model.funcs.h_obj = @prob_hobj;
model.funcs.con = @prob_con;
model.funcs.jac = @prob_jac;
model.funcs.h_con = @prob_hcon;
model.xl = [-10; -10];
model.xu = [10; 10];
model.A = sparse([1, 1]);
model.bl = 0;
model.bu = 0;
model.cl = [1.5; -10];
model.cu = [inf; inf];
% define correct solution
sol.status = 'optimal';
sol.x = [1.22474487;-1.22474487];
sol.f = 5.2768479;
end


%% callback functions for NLP problem
% objective function
function f = prob_obj(x, s)
f = exp(x(1)) * (4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + 1); 
end

% gradient of objective
function g = prob_grad(x, s)
g = exp(x(1)) * [4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+8*x(1)+6*x(2)+1;...
    4*x(1)+4*x(2)+2];
end

% hessian of objective
function H = prob_hobj(x, s)
a = 4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 16*x(1) + 10*x(2) + 9;
b = 4*x(2) + 4*x(1) + 6;
H = sparse(exp(x(1)) * [a b;b 4]);
end

% constraints
function c = prob_con(x, s)
c = [ - x(1)*x(2) + x(1) + x(2); x(1)*x(2)];
end

% Jacobian of constraints
function dc = prob_jac(x, s)
dc = sparse([-x(2)+1 x(2); -x(1)+1 x(1)]');
end

% Hessian of constraints
function d2c = prob_hcon(x, lam, s)
s = [-1 1]*lam; 
d2c=sparse([ 0  s;  s  0]); 
end
