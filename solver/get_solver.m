function [optim, param, n_pts] = get_solver(fct_vec, solver, var_opt, var_fix)
% Fit a model with respect to a dataset.
%
%    This function provides a common interface for different solvers.
%    Different error metrics (i.e. norms) are available.
%    The variable ranges are normalized.

% extract the solver data
optimizer = solver.optimizer;
n_cache = solver.n_cache;
tol_cache = solver.tol_cache;
use_cache = solver.use_cache;
vec_cache = solver.vec_cache;
tol_bound = solver.tol_bound;

% get the normalized initial points, lower bound, and upper bound
fprintf('get var\n')
obj_var = SolverVar(var_opt, var_fix, tol_bound);

% cache data
fprintf('get cache\n')
fct_vec_cache = @(x_scale) get_vec_cache(x_scale, obj_var, fct_vec);
obj_cache = SolverCache(fct_vec_cache, use_cache, vec_cache, n_cache, tol_cache);

% solver
fprintf('get interface\n')
obj_interface = SolverInterface(obj_var, obj_cache);

% run the solver
disp('run solvers')
x_scale = obj_var.get_x0_scale();
for i=1:length(optimizer)
    [x_scale, optim{i}] = obj_interface.get_run(x_scale, optimizer{i});
end
[n_pts, param] = obj_var.get_param(x_scale);

end

function [err_vec, wgt_vec] = get_vec_cache(x_scale, obj_var, fct_vec)

[n_pts, param] = obj_var.get_param(x_scale);
[err_vec, wgt_vec] = fct_vec(param, n_pts);

err_vec = err_vec.';
wgt_vec = wgt_vec.';
assert(size(err_vec, 1)==n_pts, 'invalid size: err_vec')
assert(size(wgt_vec, 1)==n_pts, 'invalid size: wgt_vec')

end