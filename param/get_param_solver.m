function solver = get_param_solver()
% Select the fitting option and variables.

% assign solver parameters
solver.use_cache = true;
solver.vec_cache = true;
solver.n_cache = 1e3;
solver.tol_cache = 1e-12;

solver.tol_bound = 0.05;
solver.error_norm = 8.0;

solver.optimizer = {};
solver.optimizer{end+1} = get_optimizer('init');
solver.optimizer{end+1} = get_optimizer('ga');
solver.optimizer{end+1} = get_optimizer('fminsearch');

end

function optimizer = get_optimizer(solver_type)
% Select the fitting option and variables.

% fitting options (MATLAB optimizer option)
switch solver_type
    case 'init'
        options.n_batch = 500;
        options.n_tot = 50;
        options.n_eval_max = 5e3;
        options.val_lim = 0.3;
                
        clamp_bnd = false;
        recover_val = NaN;
    case 'fminunc'
        options = optimset ('fminunc');
        options = optimset(options, 'TolFun', 1e-6);
        options = optimset(options, 'TolX', 1e-6);
        options = optimset(options, 'MaxFunEvals', 5e3);
        
        clamp_bnd = true;
        recover_val = NaN;
    case 'fminsearch'
        options = optimset('fminsearch');
        options = optimset(options, 'TolFun', 1e-6);
        options = optimset(options, 'TolX', 1e-6);
        options = optimset(options, 'MaxFunEvals', 5e3);
        
        clamp_bnd = true;
        recover_val = NaN;
    case 'fmincon'
        options = optimset ('fmincon');
        options = optimset(options, 'TolFun', 1e-6);
        options = optimset(options, 'TolCon', 1e-6);
        options = optimset(options, 'FinDiffRelStep', 1e-5);
        options = optimset(options, 'MaxFunEvals', 5e3);
        options = optimset(options, 'UseParallel', false);

        clamp_bnd = false;
        recover_val = NaN;
    case 'surrogateopt'
        options = optimoptions (@surrogateopt);
        options = optimoptions(options, 'ConstraintTolerance', 1e-6);
        options = optimoptions(options, 'MinSurrogatePoints', 50);
        options = optimoptions(options, 'MaxFunctionEvaluations', 5e3);
        options = optimoptions(options, 'MaxTime', 60.*60);
        options = optimoptions(options, 'UseVectorized', true);
        options = optimoptions(options, 'UseParallel', false);

        clamp_bnd = false;
        recover_val = 100.0;
    case 'particleswarm'
        options = optimoptions (@particleswarm);
        options = optimoptions(options, 'FunctionTolerance', 1e-6);
        options = optimoptions(options, 'MaxIterations', 50);
        options = optimoptions(options, 'MaxStallIterations', 50);
        options = optimoptions(options, 'SwarmSize', 500);
        options = optimoptions(options, 'MaxStallTime', 60.*60);
        options = optimoptions(options, 'MaxTime', 60.*60);
        options = optimoptions(options, 'Vectorized', 'on');
        options = optimoptions(options, 'UseParallel', true);

        clamp_bnd = false;
        recover_val = NaN;
    case 'ga'
        options = optimoptions (@ga);
        options = optimoptions(options, 'FunctionTolerance', 1e-6);
        options = optimoptions(options, 'ConstraintTolerance', 1e-6);
        options = optimoptions(options, 'Generations', 1000);
        options = optimoptions(options, 'MaxStallGenerations', 50);
        options = optimoptions(options, 'PopulationSize', 500);
        options = optimoptions(options, 'MaxStallTime', 60.*60);
        options = optimoptions(options, 'MaxTime', 60.*60);
        options = optimoptions(options, 'Vectorized', 'on');
        options = optimoptions(options, 'UseParallel', true);

        clamp_bnd = false;
        recover_val = NaN;
    otherwise
        error('invalid data')
end

optimizer.log_iter = true;
optimizer.log_final = true;
optimizer.clamp_bnd = clamp_bnd;
optimizer.recover_val = recover_val;
optimizer.solver_type = solver_type;
optimizer.options = options;

end
