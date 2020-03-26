function opts = mosif_opts(opts)
%% check provided optimization options and set missing ones to the defaults
% Note that non-default options for each solver should be provided in the
% solver-specific format, which the user can find by referring to its
% manual. This is done to ensure that the user can flexibly set any options
% supported by a desired solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS
%  opts                  - structure containing various general and solver-
%                          specific options for solving the optimization 
%                          problem. If it is not provided or is empty, all 
%                          fields are initialized to their defaults. If a 
%                          field is not provided by the user, it will be 
%                          initizalized to its default. The fields are 
%                          (default values are given in brackets):
%   .solver      'ipopt' - solver name. Supported and tested options are:
%                           TOMLAB: 'snopt', 'knitro', 'minos', 'npsol'
%                           OPTI:   'ipopt' , 'clp'
%                           'mosek'
%   .printlevel      []  - level of command window printing ([] or 0):
%                           [] - let each solver use its default
%                           0  - nothing is printed.
%   .ipopt           []  - structure whose fields are options for IPOPT. 
%   .clp                 - structure whose fields are options for CLP. 
%   .mosek               - structure whose fields are options for MOSEK. 
%   .tomlab              - structure whose fields are options for solvers 
%                          provided as part of Tomlab distribution. The
%                          fields are:
%     .sol           []  - row vector whose entries are options of Tomlab's 
%                          SOL solvers (snopt, minos, npsol). 
%     .knitro        []  - structure whose fields are options for Tomlab's 
%                          KNITRO.
%OUTPUTS
%  opts                  - structure of options with all missing fields
%                          initialized to their defaults.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fill out missing fields of options structure with empty or default values
if nargin == 0
    opts = struct();
elseif ~isstruct(opts)
    error('Variable with optimization options must be a structure or empty.');
end
opts = populate_options(opts);

end


% add option fields
function options = populate_options(options)
    options = populate_field(options,'solver','ipopt');
    options = populate_field(options,'printlevel',[]);
    options = populate_field(options,'SOL',[]);
    options = populate_field(options,'KNITRO',[]);
    options = populate_field(options,'IPOPT',[]);
    options = populate_field(options,'CLP',[]);
    options = populate_field(options,'MOSEK',[]);
end

% populate given field of a structure
function s = populate_field(s,fieldname, fieldvalue)
    if ~isfield(s,fieldname) || isempty(s.(fieldname))
        s.(fieldname)=fieldvalue;
    end
end


