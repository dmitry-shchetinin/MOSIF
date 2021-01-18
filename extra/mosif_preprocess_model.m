function model = mosif_preprocess_model(model, solver)


%% fill out missing fields of model structure with empty or default values
model = populate_model(model);

%% check model consistency
mosif_check_model(model, solver);

%% analyze and record optimization problem type
model.isNLP = ~isempty(model.funcs.obj);

%% standardize values of model fields (make matrices sparse, etc.)
model = mosif_standardize_model(model);

%% estimate sparsity patterns if needed
% Jacobian sparsity pattern
if model.isNLP && isempty(model.Js) && ~isempty(model.funcs.con)
    model.Js = Jacobian_sparsity_pattern(model);
end
% Hessian sparsity pattern
if model.isNLP && isempty(model.Hs) && ~isempty(model.funcs.hess_obj)
    model.Hs = Hessian_sparsity_pattern(model);
elseif ~model.isNLP && nnz(model.Q) > 0
    model.Hs = mosif_set_entries_to_ones(model.Q);
end

%% if model is NLP, prepare it for solution by NLP solver
NLP_solvers = ["ipopt", "knitro", "snopt"];
if any(NLP_solvers == solver)
   model = mosif_prepare_for_NLP_solver(model); 
end

end


%% add model fields
function model = populate_model(model)
    if ~isfield(model,'x0') || isempty(model.x0)
        error('Starting point must be provided.');
    end
    n = length(model.x0);
    model = populate_field(model, 'Q', sparse(n,n));
    model = populate_field(model, 'q', zeros(n,1));
    model = populate_field(model, 'p', 0);
    model = populate_field(model, 'A', []);
    model = populate_field(model, 'bl', []);
    model = populate_field(model, 'bu', []);
    model = populate_field(model, 'xl', []);
    model = populate_field(model, 'xu', []);
    model = populate_field(model, 'x0', []);
    model = populate_field(model, 'xtype', repelem('C',n));
    model = populate_field(model, 'cl', []);
    model = populate_field(model, 'cu', []);
    model = populate_field(model, 'Js', []);
    model = populate_field(model, 'Js_quadcon', []);
    model = populate_field(model, 'Hs', []);
    model = populate_field(model, 'funcs', []);
    model = populate_field(model, 'name', 'Optimization problem');
    model = populate_field(model, 'quadcon', []);
    model = populate_field(model, 'cones', []);
    model = populate_field(model, 'userdata', []);
    model.funcs = populate_field(model.funcs, 'obj', []);
    model.funcs = populate_field(model.funcs, 'grad', []);
    model.funcs = populate_field(model.funcs, 'hess_obj', []);
    model.funcs = populate_field(model.funcs, 'con', []);
    model.funcs = populate_field(model.funcs, 'jac', []);
    model.funcs = populate_field(model.funcs, 'hess_con', []);
end


%% populate given field of a structure
function s = populate_field(s, fieldname, fieldvalue)
    if ~isfield(s,fieldname) || isempty(s.(fieldname))
        s.(fieldname)=fieldvalue;
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




