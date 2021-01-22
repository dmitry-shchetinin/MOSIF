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
    quadcon = auxdata.quadcon;
    n_quadcon = length(quadcon);
    n = length(x);
    % get part of Jacobian related to linear part of constraints
    Jquadcon = [quadcon.qc]';
    % get part of Jacobian related to quadratic part of constraints
    Jvals = cell(n_quadcon, 1);
    for i = 1:n_quadcon
        Jvals{i} = 2 * quadcon(i).Qc * x(quadcon(i).x_idx);
    end
    % combine the two parts together
    Jquadcon = Jquadcon + sparse([quadcon.Jrows], [quadcon.x_idx], ...
                                 cell2mat(Jvals), n_quadcon, n);
else
    Jquadcon = [];
end

% combine all part of Jacobian
J = [Jnonlin; Jquadcon; Jlin];
end