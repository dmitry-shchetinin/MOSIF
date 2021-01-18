function model = mosif_standardize_model(model)

%% make sure all matrices are sparse
matrix_names = {'Q', 'A', 'Js', 'Hs'};
for i = 1:length(matrix_names)
    name = matrix_names{i};
    if ~issparse(model.(name))
        model.(name) = sparse(model.(name));
    end
end
% make sure quadratic coefficient matrix is symmetric
model.Q = (model.Q + model.Q') / 2; 

%% make sure that all coefficient vectors are column vectors
vector_names = {'q', 'bl', 'bu', 'cl', 'cu', 'xl', 'xu'};
for i = 1:length(vector_names)
    name = vector_names{i};
    model.(name) = make_column_vector(model.(name));
end

%% make sure all quadratic constraints have sparse matrices and vectors
for i = 1:length(model.quadcon)
    model.quadcon(i).Qc = sparse(model.quadcon(i).Qc + model.quadcon(i).Qc') / 2; % make sure it's symmetric
    model.quadcon(i).qc = sparse(model.quadcon(i).qc);
    model.quadcon(i).qc = make_column_vector(model.quadcon(i).qc);
end

%% make sure all conic constraints are defined by column vectors
for i = 1:length(model.cones)
    model.cones(i).idx_le = make_column_vector(model.cones(i).idx_le);
    model.cones(i).idx_ge = make_column_vector(model.cones(i).idx_ge);
end

%% set infinite lower and upper constraint and variable bounds if not given
% variable bounds
nv = length(model.x0);
model = populate_field(model,'xl', -inf(nv,1));
model = populate_field(model,'xu', inf(nv,1));

% linear constraints
nl = size(model.A, 1);
model = populate_field(model,'bl', -inf(nl,1));
model = populate_field(model,'bu', inf(nl,1));

% nonlinear constraints
if ~isempty(model.cl) || ~isempty(model.cu)
    nn = max(numel(model.cl), numel(model.cu));
elseif model.isNLP && ~isempty(model.funcs.con)
    nn = numel(feval(model.funcs.con, model.x0, model));
else
    nn = 0;
end
model = populate_field(model,'cl', -inf(nn,1));
model = populate_field(model,'cu', inf(nn,1));


end


%% ensure that vector is a column vector
function v = make_column_vector(v)
    if ~iscolumn(v)
        v = transpose(v);
    end
end


%% populate given field of a structure
function s = populate_field(s, fieldname, fieldvalue)
    if ~isfield(s,fieldname) || isempty(s.(fieldname))
        s.(fieldname)=fieldvalue;
    end
end