function model = mosif_prepare_for_NLP_solver(model)

%% convert conic constraints to quadratic
model = mosif_cone2quadratic(model);

%% make changes to constraint bounds and sparsity patterns
nn_quadcon = length(model.quadcon);
n = length(model.x0);
if nn_quadcon > 0
    % adjust bounds on nonlinear constraints
    model.cl = [model.cl; -inf(nn_quadcon, 1)];
    model.cu = [model.cu; [model.quadcon.rhs]'];
    % adjust Jacobian sparsity pattern
    Js_quadcon = abs([model.quadcon.qc]');
    for i = 1:nn_quadcon
        Js_quadcon(i,:) = Js_quadcon(i,:) + sum(abs(model.quadcon(i).Qc));
    end
    model.Js_quadcon = mosif_set_entries_to_ones(Js_quadcon);
    model.Js = [model.Js; model.Js_quadcon];
    % adjust Hessian sparsity pattern
    if isempty(model.Hs)
        Hs_quadcon = sparse(n,n);
    else
        Hs_quadcon = model.Hs;
    end
    for i = 1:nn_quadcon
        Hs_quadcon = Hs_quadcon + abs(model.quadcon(i).Qc);
    end
    model.Hs = mosif_set_entries_to_ones(Hs_quadcon);
end

end