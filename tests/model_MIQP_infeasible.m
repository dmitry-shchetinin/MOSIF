function [model, sol] = model_MIQP_infeasible()
%% define model of a simple infeasible MIQP problem
% the following problem is formulated
% min x1^2 + x2 + x3 + x4
% s.t. x3 - 2*x4 = -0.5
%      -1 <= x_i <= 1 (for i=1..4)
%      x3, x4 in {0,1}
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
model.name = 'MIQP_optimal';
model.x0 = zeros(4,1);
model.Q = sparse(1,1,1,4,4);
model.q = [0; 1; 1; 1];
model.A = sparse([0, 0, 1, -2]);
model.bl = -0.5;
model.bu = -0.5;
model.xl = -ones(4,1);
model.xu = ones(4,1);
model.xtype = 'CCBB';

%% define correct solution
sol.status = 'infeasible';
sol.x = [];
sol.f = [];
end