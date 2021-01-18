function J = mosif_jacobian(x, auxdata)
% from linear constraints
Jlin = auxdata.A;

% from nonlinear general constraints
if auxdata.isNLP && ~isempty(auxdata.funcs.jac)
    Jnonlin = feval(auxdata.funcs.jac, x, auxdata);
else
    Jnonlin = [];
end

% from quadratic constraints
if ~auxdata.isNLP && ~isempty(auxdata.quadcon)
    % in order not to change sparsity pattern each time, we'll use a trick
    % here: first create jacobian from a sparsity pattern, and then
    % subtract the sparsity pattern so that we end up with the actual
    % jacobian matrix
    Jquadcon = auxdata.Js_quadcon;
    n_quadcon = length(auxdata.quadcon);
    n = length(x);
    for i = 1:n_quadcon
        tmp = (2 * sparse(auxdata.quadcon(i).Qc * x) + auxdata.quadcon(i).qc)';
        [~, ic, v] = find(tmp);
        Jquadcon = Jquadcon + sparse(i*ones(length(v),1), ic, v, n_quadcon, n);
    end
    Jquadcon = Jquadcon - auxdata.Js_quadcon;
else
    Jquadcon = [];
end

% combine
J = [Jnonlin; Jquadcon; Jlin];
end