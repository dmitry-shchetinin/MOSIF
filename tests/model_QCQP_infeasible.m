function [model, sol] = model_QCQP_infeasible()
%% define model of a simple infeasible QCQP problem
% the following problem is formulated
% min x1^2 + x2 + x3 + x4
% s.t. x3 - 2*x4 = 0
%      -5 <= x_i <= 3 (for i=1..4)
%      x2^2 - 2*x2 <= 0
%      x4^2 - 4*x4 <= -4
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
model.name = 'QCQP_infeasible';
model.x0 = zeros(4,1);
model.Q = sparse(1,1,1,4,4);
model.q = [0; 1; 1; 1];
model.A = sparse([0, 0, 1, -2]);
model.bl = 0;
model.bu = 0;
model.xl = -5 * ones(4,1);
model.xu = 3 * ones(4,1);
model.quadcon(1).Qc = sparse(2,2,1,4,4);
model.quadcon(1).qc = sparse(2,1,-2,4,1);
model.quadcon(1).rhs = 0;
model.quadcon(2).Qc = sparse(4,4,1,4,4);
model.quadcon(2).qc = sparse(4,1,-4,4,1);
model.quadcon(2).rhs = -4;

%% define correct solution
sol.status = 'infeasible';
sol.x = [];
sol.f = [];
end