function c = mosif_constraints(x, auxdata)
% add linear constraints
if ~isempty(auxdata.A)
    clin = auxdata.A * x;
else
    clin = [];
end

% add nonlinear general constraints
if auxdata.isNLP && ~isempty(auxdata.funcs.con)
    cnonlin = feval(auxdata.funcs.con, x, auxdata);
else
    cnonlin = [];
end

% add quadratic constraints
if ~auxdata.isNLP && ~isempty(auxdata.quadcon)
    quadcon = auxdata.quadcon;
    n_quadcon = length(quadcon);
    cquadcon = zeros(n_quadcon, 1);
    for i = 1:n_quadcon
        xc = x(quadcon(i).x_idx); % entries of x present in quadratic part of the constraint
        cquadcon(i) = xc' * quadcon(i).Qc * xc + quadcon(i).qc' * x;
    end
else
    cquadcon = [];
end

% combine into one vector
c = [cnonlin; cquadcon; clin];
end