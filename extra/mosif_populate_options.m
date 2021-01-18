function opts = mosif_populate_options(opts)
%% check provided optimization options and set missing ones to the defaults
% Note that non-default options for each solver should be provided in the
% solver-specific format, which the user can find by referring to its
% manual. This is done to ensure that the user can flexibly set any options
% supported by a desired solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%  opts               - structure containing various general and solver-
%                       specific options for solving the optimization 
%                       problem. For details, see function 'solve_optprob'.
%OUTPUTS
%  opts               - structure of options with all missing fields
%                       initialized to their defaults.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check that options is a structure
if ~isstruct(opts)
    error('Variable with optimization options must be a structure or empty.');
end

%% fill out missing fields of the structure with empty or default values
opts = populate_field(opts, 'solver', 'ipopt');
opts = populate_field(opts, 'suppress_print', 0);
opts = populate_field(opts, 'ipopt', []);
opts = populate_field(opts, 'clp', []);
opts = populate_field(opts, 'gurobi', []);
opts = populate_field(opts, 'mosek', []);
opts = populate_field(opts, 'tomlab', []);
opts.tomlab = populate_field(opts.tomlab, 'knitro', []);
opts.tomlab = populate_field(opts.tomlab, 'sol', []);
end


%% populate given field of a structure
function s = populate_field(s, fieldname, fieldvalue)
    if ~isfield(s,fieldname) || isempty(s.(fieldname))
        s.(fieldname)=fieldvalue;
    end
end


