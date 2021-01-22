function model = mosif_prepare_for_NLP_solver(model)

%% estimate sparsity patterns of derivatives of nonlinear constraints
if model.isNLP
    % Jacobian sparsity pattern
    if isempty(model.Js) && ~isempty(model.funcs.con)
        model.Js = Jacobian_sparsity_pattern(model);
    end
    % Hessian sparsity pattern
    if isempty(model.Hs) && ~isempty(model.funcs.hess_obj)
        model.Hs = Hessian_sparsity_pattern(model);
    end
    return;
end

%% get sparsity pattern of quadratic objective
if nnz(model.Q) > 0
    model.Hs = mosif_set_entries_to_ones(model.Q);
end

%% convert conic constraints to quadratic
model = mosif_cone2quadratic(model);

%% if there are no quadratic constraints, there is nothing to adjust
n_quadcon = length(model.quadcon);
if n_quadcon == 0
    return;
end
n = length(model.x0);

%% adjust bounds on nonlinear constraints
model.cl = [model.cl; -inf(n_quadcon, 1)];
model.cu = [model.cu; [model.quadcon.rhs]'];

%% adjust representation of nonlinear constraints
% this is done to enable faster evaluation of constraints and their
% derivatives in callback functions of NLP solver
for i = 1:n_quadcon
    % get nonzeros indices and values of Qc matrix
    [Hrows, Hcols, Hvals] = find(model.quadcon(i).Qc);
    model.quadcon(i).Hrows = Hrows';
    model.quadcon(i).Hcols = Hcols';
    model.quadcon(i).Hvals = Hvals';
    % get indices of entries of vector x that correspond to QC nonzeros
    x_idx = unique(Hrows);
    model.quadcon(i).x_idx = x_idx';
    % reduce matrix Qc to its nonzero block
    model.quadcon(i).Qc = sparse(model.quadcon(i).Qc(x_idx, x_idx));
    % record row indices for part of sparse matrix that defines i-th
    % constraint in Jacobian of quadratic constraints
    model.quadcon(i).Jrows = i * ones(1, length(x_idx));
end
    
%% get Jacobian sparsity pattern
% part related to linear part of quadratic constraints
model.Js = abs([model.quadcon.qc]');
% add part related to quadratic part of linear constraints
rows = [model.quadcon.Jrows];
cols = [model.quadcon.x_idx];
vals = ones(length(rows), 1);
model.Js = model.Js + sparse(rows, cols, vals, n_quadcon, n);
model.Js = mosif_set_entries_to_ones(model.Js);

%% get Hessian sparsity pattern
rows = [model.quadcon.Hrows];
cols = [model.quadcon.Hcols];
vals = ones(length(rows), 1);
Hs_quadcon = sparse(rows, cols, vals, n, n);
if isempty(model.Hs)
    model.Hs = mosif_set_entries_to_ones(Hs_quadcon);
else
    model.Hs = mosif_set_entries_to_ones(model.Hs + Hs_quadcon);
end
  

end



%% estimate Jacobian sparsity pattern
function Js = Jacobian_sparsity_pattern(model)
    % get number of variables
    n = numel(model.x0);
    % compute sum of absolute values of Jacobian for 3 random values of x
    Js = abs(feval(model.funcs.jac, rand(n,1), model)) + ...
        abs(feval(model.funcs.jac, rand(n,1), model)) + ...
        abs(feval(model.funcs.jac, rand(n,1), model));
    Js = mosif_set_entries_to_ones(sparse(Js));
end


%% estimate Hessian sparsity pattern
function Hs = Hessian_sparsity_pattern(model)
    % get number of variables
    nv = numel(model.x0);
    % compute sum of absolute values of Hessian for 3 random values of x
    Hs = abs(feval(model.funcs.hess_obj, rand(nv,1), model))+...
        abs(feval(model.funcs.hess_obj, rand(nv,1), model))+...
        abs(feval(model.funcs.hess_obj, rand(nv,1), model));
    if ~isempty(model.funcs.hess_con)
        % get number of nonlinear constraints
        nn = numel(feval(model.funcs.con, model.x0, model));
        % compute sum of absolute values of Hessian for 3 random values of x
        Hs=Hs+abs(feval(model.funcs.hess_con, rand(nv,1), rand(nn,1), model)) + ...
            abs(feval(model.funcs.hess_con, rand(nv,1), rand(nn,1), model)) + ...
            abs(feval(model.funcs.hess_con, rand(nv,1), rand(nn,1), model));
    end
    Hs = mosif_set_entries_to_ones(sparse(Hs));
end
