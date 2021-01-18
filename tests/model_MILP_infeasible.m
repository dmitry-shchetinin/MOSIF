function [model, sol] = model_MILP_infeasible()
%% define model of a simple infeasible MILP problem
% the following problem is formulated
% min  x1 + 2*x2 + 4*x3 + 8*x4
% s.t. x1 + x2 + x3 + x4 = 1.5
%      0 <= x_i <= 1 (for i=1..4)
%      x1, x2, x3, x4 in {0,1}
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
model.name = 'MILP_infeasible';
model.x0 = zeros(4,1);
model.q = [1; 2; 4; 8];
model.xl = zeros(4,1);
model.xu = ones(4,1);
model.xtype = 'BBBB';
model.A = sparse(ones(1,4));
model.bl = 1.5;
model.bu = 1.5;

%% define correct solution
sol.status = 'infeasible';
sol.x = [];
sol.f = [];
end