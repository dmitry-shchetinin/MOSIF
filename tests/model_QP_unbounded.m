function [model, sol] = model_QP_unbounded()   % NOTE: clp fails on this
%% define model of a simple unbounded QP problem
% the following problem is formulated
% min x1^2 + x2 + x3 + x4
% s.t. x1 = 1 (as a constraint)
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
model.name = 'QP_unbounded';
model.x0 = zeros(4,1);
model.Q = sparse(1,1,1,4,4);
model.q = [0; 1; 1; 1];
model.A = sparse([1, 0, 0, 0]);
model.bl = 1;
model.bu = 1;

%% define correct solution
sol.status = 'unbounded';
sol.x = [];
sol.f = [];
end