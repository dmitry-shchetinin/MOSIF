function model = mosif_cone2quadratic(model)
n_cones = length(model.cones);
n = length(model.x0);
n_quad = length(model.quadcon);
for i = 1:n_cones
    if length(model.cones(i).idx_ge) == 1
        % if it is a normal cone
        idx = [model.cones(i).idx_ge; model.cones(i).idx_le];
        model.quadcon(n_quad+i).Qc = sparse(idx, idx, [-1; ones(length(idx)-1,1)], n, n);
    else
        % if it is a rotated cone
        idx = model.cones(i).idx_le;
        idx1 = model.cones(i).idx_ge(1);
        idx2 = model.cones(i).idx_ge(2);
        model.quadcon(n_quad+i).Qc = sparse([idx1; idx2; idx], [idx2; idx1; idx], ...
                                            [-0.5; -0.5; ones(length(idx),1)], n, n);
    end
    model.quadcon(n_quad+i).qc = sparse(n,1);
    model.quadcon(n_quad+i).rhs = 0;
end
model.cones = [];
end
