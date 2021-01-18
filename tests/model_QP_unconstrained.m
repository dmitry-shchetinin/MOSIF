function [model, sol] = model_QP_unconstrained()
%% define model of a simple unconstrained QP problem
% the following problem is formulated
% min x1^2 + x2^2 + x3^2 + x4^2
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
model.name = 'QP_unconstrained';
model.x0 = [1; 2; 3; 4];
model.Q = speye(4);

%% define correct solution
sol.status = 'optimal';
sol.x = zeros(4,1);
sol.f = 0;
end