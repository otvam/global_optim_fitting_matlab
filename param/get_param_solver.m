function [cache, optimizer] = get_param_solver()
% Select the solver type and options.
%
%    Returns:
%        cache (struct): structure with the cache options
%        optimizer (cell): cell describing the different solvers to be used
%
%    Thomas Guillod.
%    2021-2022 - BSD License.

% use (or not) the cache
cache.use_cache = true;

% allow (or not) vectorized call to the error function
%    - if true, several parameters combinations are evaluated at once (vectorized call)
%    - if false, each combination is evaluated separately (in a for-loop)
cache.vec_cache = true;

% use (or not) parfor for evaluating many combinations in a for-loop
cache.parfor_cache = false;

% maximum number of elements in the cache
cache.n_cache = 1e3;

% tolerance for determining unicity in the cache
cache.tol_cache = 1e-12;

% solver to be used for solving the problem
%    - the different solvers are called sequentially
%    - the first solver is using the provided initial values
%    - afterwards, the results is used as initial values
optimizer = {};
optimizer{end+1} = get_optimizer('init');
optimizer{end+1} = get_optimizer('ga');
optimizer{end+1} = get_optimizer('fminsearch');

end

function optimizer = get_optimizer(solver_type)
% Select the solver options for a specific solver.
%
%    Parameters:
%        solver_type (str): solver name
%
%    Returns:
%        optimizer (struct): structure describing the solver parameters

switch solver_type
    case 'init'
        % brute force solver for generating reasonable initial values
        %    - generate random parameters combinations (within bounds)
        %    - select the combinations with an error metric below the threshold
        options.err_lim = 0.3; % threshold for keeping a parameter combinations
        options.n_tot = 50; % number of combinations to be found
        options.n_batch = 50; % number of tested combinations per solver iteration
        options.n_iter_max = 25; % maximum number of iterations
        
        % solver handles constraints, sine transformation is not required
        clamp_bnd = false;
        
        % solver can deal with undefined values, keep such values
        recover_val = NaN;
    case 'fminunc'
        % MATLAB fminunc options
        options = optimset ('fminunc');
        options = optimset(options, 'TolFun', 1e-6);
        options = optimset(options, 'TolX', 1e-6);
        options = optimset(options, 'MaxFunEvals', 5e3);
        
        % solver does not provide bounds, add constraints with sine transformation
        clamp_bnd = true;
        
        % solver cannot deal with undefined values, replace such values
        recover_val = 100.0;
    case 'fminsearch'
        % MATLAB fminsearch options
        options = optimset('fminsearch');
        options = optimset(options, 'TolFun', 1e-6);
        options = optimset(options, 'TolX', 1e-6);
        options = optimset(options, 'MaxFunEvals', 5e3);
        
        % solver does not provide bounds, add constraints with sine transformation
        clamp_bnd = true;
        
        % solver cannot deal with undefined values, replace such values
        recover_val = 100.0;
    case 'fmincon'
        % MATLAB fmincon options
        options = optimset ('fmincon');
        options = optimset(options, 'TolFun', 1e-6);
        options = optimset(options, 'TolCon', 1e-6);
        options = optimset(options, 'FinDiffRelStep', 1e-5);
        options = optimset(options, 'MaxFunEvals', 5e3);
        options = optimset(options, 'UseParallel', false);
        
        % solver handles constraints, sine transformation is not required
        clamp_bnd = false;
        
        % solver cannot deal with undefined values, replace such values
        recover_val = 100.0;
    case 'surrogateopt'
        % MATLAB surrogateopt options
        options = optimoptions (@surrogateopt);
        options = optimoptions(options, 'ConstraintTolerance', 1e-6);
        options = optimoptions(options, 'MinSurrogatePoints', 50);
        options = optimoptions(options, 'MaxFunctionEvaluations', 5e3);
        options = optimoptions(options, 'MaxTime', 60.*60);
        options = optimoptions(options, 'UseVectorized', true);
        options = optimoptions(options, 'UseParallel', false);
        
        % solver handles constraints, sine transformation is not required
        clamp_bnd = false;
        
        % solver cannot deal with undefined values, replace such values
        recover_val = 100.0;
    case 'particleswarm'
        % MATLAB particleswarm options
        options = optimoptions (@particleswarm);
        options = optimoptions(options, 'FunctionTolerance', 1e-6);
        options = optimoptions(options, 'MaxIterations', 50);
        options = optimoptions(options, 'MaxStallIterations', 50);
        options = optimoptions(options, 'SwarmSize', 500);
        options = optimoptions(options, 'MaxStallTime', 60.*60);
        options = optimoptions(options, 'MaxTime', 60.*60);
        options = optimoptions(options, 'Vectorized', 'on');
        options = optimoptions(options, 'UseParallel', true);
        
        % solver handles constraints, sine transformation is not required
        clamp_bnd = false;
        
        % solver can deal with undefined values, keep such values
        recover_val = NaN;
    case 'ga'
        % MATLAB ga options
        options = optimoptions (@ga);
        options = optimoptions(options, 'FunctionTolerance', 1e-6);
        options = optimoptions(options, 'ConstraintTolerance', 1e-6);
        options = optimoptions(options, 'Generations', 100);
        options = optimoptions(options, 'MaxStallGenerations', 50);
        options = optimoptions(options, 'PopulationSize', 500);
        options = optimoptions(options, 'MaxStallTime', 60.*60);
        options = optimoptions(options, 'MaxTime', 60.*60);
        options = optimoptions(options, 'Vectorized', 'on');
        options = optimoptions(options, 'UseParallel', true);
        
        % solver handles constraints, sine transformation is not required
        clamp_bnd = false;
        
        % solver can deal with undefined values, keep such values
        recover_val = NaN;
    otherwise
        error('invalid data')
end

% solver name
optimizer.solver_type = solver_type;

% log the solver progress after each iteration
optimizer.log_iter = true;

% log the solver results after the final iteration
optimizer.log_final = true;

% handling of variable bounds
%    - if false, the solver is handling the constraints
%    - if true, the contraints are implemented with a sine transformation
%    - this allows to implements bounds with solver that do not provide constraints
optimizer.clamp_bnd = clamp_bnd;

% handling of non-finite values (NaN and Inf) returned by the error function
%    - some solvers are not supporting non-finite values for the error functions
%    - the non-finite values can be replaced by a (large) constant value
optimizer.recover_val = recover_val;

% solver options
optimizer.options = options;

end
