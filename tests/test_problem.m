clear;

n = 5;
model.type = 'QP';
model.Q = speye(n);
model.c = (1:n)';
model.A = sparse(ones(1,n));
model.bl = 1;
model.bu = 1;
model.x0 = rand(n,1);
opt = [];
%model.xl = zeros(n,1);
opt.solver = 'mosek';

tic;
[result, extra] = solve_problem(model, opt);
time = toc;