function [model, sol] = model_NLP_optimal()
%% define model of a simple NLP problem that has optimal solution
% the following problem is formulated
% min exp(x1) * (4*x1^2 + 2*x2^2 + 4*x1*x2 + 2*x2 + 1)
% s.t. -10<=x1<=10
%      -10<=x2<=10
%      x1+x2=0
%      -x1*x2 + x1 + x2>=1.5
%      x1*x2>=-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUTS
%  model                - structure containing formulation of the problem
%                         in mosif format (see 'solve_optprob' for details)
%  sol                  - structure representing solution of the problem,
%                         with the following fields:
%   .f                  - optimal value of objective function
%   .x                  - optimal value of variable vector
%   .status             - character array that contains status of solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define problem
model.name = 'NLP_optimal';
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

%% define correct solution
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
