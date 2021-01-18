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
    n_quadcon = length(auxdata.quadcon);
    cquadcon = zeros(n_quadcon, 1);
    for i = 1:n_quadcon
        cquadcon(i) = x' * auxdata.quadcon(i).Qc * x + auxdata.quadcon(i).qc' * x;
    end
else
    cquadcon = [];
end

% combine into one vector
c = [cnonlin; cquadcon; clin];
end