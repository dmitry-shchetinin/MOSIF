function [model, sol] = model_LP_unbounded()
%% define model of a simple unbounded LP problem
% the following problem is formulated
% min x1 + x2 + x3 + x4
% s.t. x1 + x2 + x3 + x4 <= 0
%      x_i <= 1 (for i=1..4)
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
model.name = 'LP_unbounded';
model.x0 = zeros(4,1);
model.q = ones(4,1);
model.xu = ones(4,1);
model.A = sparse(ones(1,4));
model.bu = 0;

%% define correct solution
sol.status = 'unbounded';
sol.x = [];
sol.f = [];
end