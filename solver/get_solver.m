function [param, optim] = get_solver(fct_err, var_opt, var_fix, solver)
% Fit parameters with respect to a dataset.
%
%    This function provides a common interface for different solvers
%        - fminunc / fminsearch / fmincon
%        - surrogateopt / particleswarm / ga
%        - combination of several solvers
%
%    Customized error function
%        - custom weights for the dataset points
%        - choice of the error norm
%        - recover from undefined values
%        - vectorized evaluation of the error function
%        - caching of the error function
%
%    Advanced parameters handling
%        - scalar or vector parameters
%        - initial values
%        - bounds (even for unconstrained solvers)
%        - variable transformation
%        - variable normlization
%
%    Advanced monitoring capabilities
%        - error metrics
%        - solver metrics
%        - computational cost / timing
%        - plots
%
%    Parameters:
%        fct_err (handle): error function for determining the parameters
%        var_opt (cell): description of the parameters to be fitted
%        var_fix (cell): description of the parameters with fixed values
%        solver (struct): structure describing the solver parameters
%
%    Returns:
%        param (struct): fitted parameters
%        optim (cell): solver metrics
%
%    Thomas Guillod.
%    2021-2022 - BSD License.

% extract the solver data
optimizer = solver.optimizer;
n_cache = solver.n_cache;
tol_cache = solver.tol_cache;
use_cache = solver.use_cache;
vec_cache = solver.vec_cache;
tol_bound = solver.tol_bound;

% objected managing the variables (name, bounds, transformation, normalization, etc.)
fprintf('get var\n')
obj_var = SolverVar(var_opt, var_fix, tol_bound);

% cache object for the error function
fprintf('get cache\n')
fct_err_cache = @(x_scale) get_vec_cache(x_scale, obj_var, fct_err);
obj_cache = SolverCache(fct_err_cache, use_cache, vec_cache, n_cache, tol_cache);

% object interfacing the different solvers
fprintf('get interface\n')
obj_interface = SolverInterface(obj_var, obj_cache);

% calling the different solvers (using the results as initial values)
disp('run solvers')
x_scale = obj_var.get_x0_scale();
for i=1:length(optimizer)
    [x_scale, optim{i}] = obj_interface.get_run(x_scale, optimizer{i});
end

% extract the parameter structure from a raw matrix (transformation and normalization)
[n_pts, param] = obj_var.get_param(x_scale);

% solution should be a single parameter combination
assert(n_pts==1, 'invalid solution')

end

function [err_mat, wgt_mat] = get_vec_cache(x_scale, obj_var, fct_err)
% Error function used by the solver (throught the cache).

% extract the parameter structure from a raw matrix (transformation and normalization)
[n_pts, param] = obj_var.get_param(x_scale);

% call the provided error function
[err_mat, wgt_mat] = fct_err(param, n_pts);

% reshape and check size
err_mat = err_mat.';
wgt_mat = wgt_mat.';
assert(size(err_mat, 1)==n_pts, 'invalid size: err_mat')
assert(size(wgt_mat, 1)==n_pts, 'invalid size: wgt_mat')

end