function model = mosif_preprocess_model(model, solver)

%% create list of supported solvers
solvers.nlp = ["ipopt", "knitro", "snopt"];
solvers.qcqp = ["gurobi", "mosek", "cplex", "ecos"];
solvers.mip = ["gurobi", "mosek", "cplex"];

%% fill out missing fields of model structure with empty or default values
model = populate_model(model);

%% check model consistency
mosif_check_model(model, solvers, solver);

%% analyze and record optimization problem type
model.isNLP = ~isempty(model.funcs.obj);

%% standardize values of model fields (make matrices sparse, etc.)
model = mosif_standardize_model(model);

%% if model is NLP, prepare it for solution by NLP solver
if any(solvers.nlp == solver)
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






