function [model, sol] = model_LP_infeasible()
%% define model of a simple infeasible LP problem
% the following problem is formulated
% min x1 + x2 + x3 + x4
% s.t. x1 + x2 + x3 + x4 = 1
%      x_i <= 0 (for i=1..4)
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
model.name = 'LP_infeasible';
model.x0 = zeros(4,1);
model.q = ones(4,1);
model.xu = zeros(4,1);
model.A = sparse(ones(1,4));
model.bl = 1;
model.bu = 1;

%% define correct solution
sol.status = 'infeasible';
sol.x = [];
sol.f = [];
end