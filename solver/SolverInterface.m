classdef SolverInterface < handle
    %% properties
    properties (SetAccess = private, GetAccess = private)
        obj_var
        obj_cache
    end
    
    %% public
    methods (Access = public)
        function self = SolverInterface(obj_var, obj_cache)
            % set data
            self.obj_var = obj_var;
            self.obj_cache = obj_cache;
            
            % disable warning
            warning('off', 'optimlib:commonMsgs:NoPCTLicense');
        end
        
        function [x_scale, optim] = get_run(self, x0_scale, optimizer)
            % Common interface for different solvers.
                        
            % extract
            solver_type = optimizer.solver_type;
            clamp_bnd = optimizer.clamp_bnd;
            recover_val = optimizer.recover_val;
            error_norm = optimizer.error_norm;
            options = optimizer.options;
            
            % scale
            self.obj_cache.get_clear();
            fct_unclamp = @(x0_scale) self.obj_var.get_unclamp(x0_scale, clamp_bnd);
            fct_clamp = @(x_clamp) self.obj_var.get_clamp(x_clamp, clamp_bnd);
            fct_param = @(x_scale) self.obj_var.get_param(x_scale);
            fct_opt = @(x_scale) self.obj_cache.get_eval(x_scale);
            
            % objective
            fct_sol = @(x_unclamp) SolverInterface.get_sol(x_unclamp, fct_opt, fct_clamp, error_norm, recover_val);
            
            % get logging
            data_optim.fct_opt = fct_opt;
            data_optim.fct_clamp = fct_clamp;
            data_optim.fct_param = fct_param;
            data_optim.error_norm = error_norm;
            data_optim.n_var = size(x0_scale, 2);
            
            % call the solver
            [x_scale, optim] = SolverInterface.get_solver(fct_sol, fct_unclamp, fct_clamp, x0_scale, options, solver_type, data_optim);
        end
    end
    
    %% private api
    methods (Static, Access = private)
        function [x_scale, optim] = get_solver(fct_sol, fct_unclamp, fct_clamp, x0_scale, options, solver_type, data_optim)
            % Call the solver.
            
            % get logger
            obj_log = SolverLog(solver_type, data_optim);
            fct_iter = @(x_unclamp, n_iter, n_eval, msg) obj_log.get_iter(x_unclamp, n_iter, n_eval, msg);
            fct_final = @(x_unclamp, n_iter, n_eval, msg, is_valid) obj_log.get_final(x_unclamp, n_iter, n_eval, msg, is_valid);
            
            % umclamp data
            [x0_unclamp, lb_unclamp, ub_unclamp] = fct_unclamp(x0_scale);
            
            % run the solver
            x_unclamp = SolverList.get_solver(fct_sol, fct_iter, fct_final, x0_unclamp, lb_unclamp, ub_unclamp, options, solver_type);
            
            % clamp data
            x_scale = fct_clamp(x_unclamp);
            
            % get solution struct
            optim = obj_log.get_optim();
        end
        
        function err = get_sol(x_unclamp, fct_opt, fct_clamp, error_norm, recover_val)
            % Objective function.
            
            % evaluate points
            x_scale = fct_clamp(x_unclamp);
            [err_vec, wgt_vec] = fct_opt(x_scale);
                                                                        
            % get error
            err = SolverUtils.get_norm(err_vec, wgt_vec, error_norm);
                        
            % recover from bad values
            idx = isfinite(err)==false;
            err(idx) = recover_val;
        end
    end
end