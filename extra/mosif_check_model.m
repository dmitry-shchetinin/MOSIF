function model = mosif_check_model(model, solvers, solver)

%% check if the solver exists
if exist(solver, 'file') ~= 3
    error('Solver %s does not exist or is not in the path.', solver)
end

%% check if the solver is suitable for the given problem
% create lists of different solvers
if ~isempty(model.funcs.obj) && all(solvers.nlp ~= solver)
    error('Solver %s cannot handle NLP problems', solver);
elseif ~isempty(model.quadcon) && all(solvers.qcqp ~= solver) && all(solvers.nlp ~= solver)
    error('Solver %s cannot handle quadratically-constrained problems', solver);
elseif any(model.xtype ~= 'C') && all(solvers.mip ~= solver)
    error('Solver %s cannot handle MIP problems', solver);
end


%% check the consistency of linear constraints
if ~isempty(model.A) && numel(model.x0) ~= size(model.A, 2)
    error('Number of columns in the matrix of linear constrains differs from number of variables');
end
 

% and so on


end

