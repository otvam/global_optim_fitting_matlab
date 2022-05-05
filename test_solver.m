function test_solver()

close('all')
addpath('solver')

solver = get_param_solver();
[var_opt, var_fix, fct_vec] = get_param_problem();

[optim, param] = get_solver(fct_vec, solver, var_opt, var_fix);

end

