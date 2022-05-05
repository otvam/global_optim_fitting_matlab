function test_solver()

close('all')
addpath('solver')

solver = get_param_solver();
[var_opt, var_fix, fct_err] = get_param_problem();

[param, optim] = get_solver(fct_err, var_opt, var_fix, solver);

end

