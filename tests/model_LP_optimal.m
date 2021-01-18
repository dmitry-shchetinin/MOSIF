function [model, sol] = model_LP_optimal()
%% define model of a simple LP problem that has optimal solution
% the following problem is formulated
% min x1 + 2*x2 + 3*x3 + 4*x4
% s.t. x1 + x2 + x3 + x4 = 1
%      0 <= x_i (for i=1..4)
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
model.name = 'LP_optimal';
model.x0 = zeros(4,1);
model.q = [1; 2; 3; 4];
model.xl = zeros(4,1);
model.A = sparse(ones(1,4));
model.bl = 1;
model.bu = 1;

%% define correct solution
sol.status = 'optimal';
sol.x = [1; 0; 0; 0];
sol.f = 1;
end