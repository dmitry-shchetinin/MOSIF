function H = mosif_hessian(x, sigma, lambda, auxdata)
H = sparse(length(x), length(x));
% add hessian of objective
if ~auxdata.isNLP
    H = H + 2 * sigma * auxdata.Q;
elseif ~isempty(auxdata.funcs.hess_obj)
    H = H + sigma * feval(auxdata.funcs.hess_obj, x, auxdata);
end

% add hessian of constraints
if ~isempty(auxdata.A)
    nl = size(auxdata.A, 1); % number of linear constraints
    lambda = lambda(1:end - nl);
end
if ~auxdata.isNLP && ~isempty(auxdata.quadcon)
    n_quadcon = length(auxdata.quadcon);
    for i = 1:n_quadcon
        H = H + 2 * lambda(i) * auxdata.quadcon(i).Qc;
    end
elseif auxdata.isNLP && ~isempty(auxdata.funcs.hess_con)
    H = H + feval(auxdata.funcs.hess_con, x, lambda, auxdata);
end

% get lower-triangular part
H = tril(H);
end