classdef SolverLog < handle
    % Class for displaying and logging the solver progress / results.
    %
    %    Log and display the solver progress after each iteration.
    %    Log, display, and plot the solver results after the final iteration.
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.

    %% properties
    properties (SetAccess = private, GetAccess = private)
        solver_type % name of the used solver
        data_optim % data required for the logging (error function and variable scaling)
        
        t_start % timestamp set at the initialization
        t_last % timestamp of the last iteration
        optim % cell containing the logged data
    end
    
    %% public
    methods (Access = public)
        function self = SolverLog(solver_type, data_optim)
            % Constructor.
            
            % set data
            self.data_optim = data_optim;
            self.solver_type = solver_type;
            
            % init the timing and logging data
            self.t_start = datetime('now');
            self.t_last = datetime('now');
            self.optim = {};
        end
        
        function get_iter(self, x_unclamp, err, n_iter, n_eval, msg)
            % Log and display the solver progress after an iteration.
            
            % get the total solver time and the time of the last iteration
            t_now = datetime('now');
            t_solver = t_now-self.t_start;
            t_iter = t_now-self.t_last;
            self.t_last = t_now;
            
            % get message
            is_valid = true;
            msg = sprintf('intermediate results\nstate: %s', msg);
            
            % get and add the logging data
            optim_tmp = SolverLog.get_log(x_unclamp, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter, self.data_optim);
            self.optim{end+1} = optim_tmp;
            
            % display the data
            name = sprintf('iter / %s / n = %d / %d', self.solver_type, n_iter, n_eval);
            SolverLog.get_disp(name, self.optim{end})
        end
        
        function get_final(self, x_unclamp, err, n_iter, n_eval, msg, is_valid)
            % Log, display, and plot the solver results after the final iteration.
            
            % get the total solver time and the average time per iteration
            t_now = datetime('now');
            t_solver = t_now-self.t_start;
            t_iter = t_solver./(n_iter+1);
            
            % get message
            msg = sprintf('final results\n%s', msg);
            
            % get and add the logging data
            optim_tmp = SolverLog.get_log(x_unclamp, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter, self.data_optim);
            self.optim{end+1} = optim_tmp;
            
            % display and plot the data
            name = sprintf('final / %s', self.solver_type);
            SolverLog.get_disp(name, self.optim{end})
            SolverLog.get_plot_single(name, self.optim{end})
            SolverLog.get_plot_all(name, self.optim)
        end
        
        function optim = get_optim(self)
            % Get the logged data.
            
            optim = self.optim;
        end
    end
    
    %% public static api
    methods(Static, Access = public)
        function get_disp(name, optim)
            % Display the logged data for a specific iteration.
            
            % name of the iteration
            fprintf('    %s\n', name)
            
            % displaty the solver message
            fprintf('        msg\n')
            for i=1:length(optim.msg)
                fprintf('            %s\n', optim.msg{i})
            end
            
            % displat the solver figures of merit
            fprintf('        fom\n')
            fprintf('            status\n')
            fprintf('                is_valid = %s\n',  mat2str(optim.sol_fom.is_valid))
            fprintf('                is_bound = %s\n',  mat2str(optim.sol_fom.is_bound))
            fprintf('            count\n')
            fprintf('                n_iter = %d\n', optim.sol_fom.n_iter)
            fprintf('                n_eval = %d\n', optim.sol_fom.n_eval)
            fprintf('            timing\n')
            fprintf('                t_solver = %s\n', char(optim.sol_fom.t_solver))
            fprintf('                t_iter = %s\n', char(optim.sol_fom.t_iter))
            fprintf('            error\n')
            fprintf('                err = %.3f %%\n', 1e2.*optim.sol_fom.err)
            fprintf('                n_pop_all = %d\n', optim.sol_fom.n_pop_all)
            fprintf('                n_pop_fail = %d\n', optim.sol_fom.n_pop_fail)
            
            % display the error metrics
            fprintf('        err_fom\n')
            fprintf('            size\n')
            fprintf('                n_rep = %d\n', optim.err_fom.n_rep)
            fprintf('                n_all = %d\n', optim.err_fom.n_all)
            fprintf('            norm\n')
            fprintf('                norm_avg = %.3f %%\n', 1e2.*optim.err_fom.norm_avg)
            fprintf('                norm_n_2 = %.3f %%\n', 1e2.*optim.err_fom.norm_n_2)
            fprintf('                norm_n_4 = %.3f %%\n', 1e2.*optim.err_fom.norm_n_4)
            fprintf('                norm_n_6 = %.3f %%\n', 1e2.*optim.err_fom.norm_n_6)
            fprintf('                norm_n_8 = %.3f %%\n', 1e2.*optim.err_fom.norm_n_8)
            fprintf('                norm_n_10 = %.3f %%\n', 1e2.*optim.err_fom.norm_n_10)
            fprintf('                norm_n_12 = %.3f %%\n', 1e2.*optim.err_fom.norm_n_12)
            fprintf('                norm_inf = %.3f %%\n', 1e2.*optim.err_fom.norm_inf)
            fprintf('            percentile\n')
            fprintf('                percentile_50 = %.3f %%\n', 1e2.*optim.err_fom.percentile_50)
            fprintf('                percentile_75 = %.3f %%\n', 1e2.*optim.err_fom.percentile_75)
            fprintf('                percentile_90 = %.3f %%\n', 1e2.*optim.err_fom.percentile_90)
            fprintf('                percentile_95 = %.3f %%\n', 1e2.*optim.err_fom.percentile_95)
            fprintf('                percentile_99 = %.3f %%\n', 1e2.*optim.err_fom.percentile_99)
            
            % display the current best parameter combinations
            fprintf('        param\n')
            field = fieldnames(optim.param);
            for i=1:length(field)
                value = optim.param.(field{i});
                fprintf('            param.%s = %s;\n', field{i}, SolverUtils.get_disp_vec(value))
            end
            
            % display if the parameters are close to the bounds
            fprintf('        bnd\n')
            field = fieldnames(optim.bnd);
            for i=1:length(field)
                value = optim.bnd.(field{i});
                fprintf('            bnd.%s = %s;\n', field{i}, SolverUtils.get_disp_vec(value))
            end
        end
        
        function get_plot_single(name, optim)
            % Plot the logged data for a specific iteration.
            
            % extract the error vector and the weighted error vector
            err_wgt_vec = optim.err_wgt_vec;
            err_vec = optim.err_vec;
                                    
            % plot the error distribution
            figure('name', [name ' / single'])
            histogram(1e2.*err_wgt_vec, 'Normalization', 'pdf')
            hold('on')
            histogram(1e2.*err_vec, 'Normalization', 'pdf')
            grid('on')
            xlabel('err (%)')
            ylabel('p (#)')
            legend('weighted', 'raw')
            title(sprintf('Histogram / %s / n = %d / %d', name, length(err_wgt_vec), length(err_vec)), 'interpreter', 'none')
        end
        
        function get_plot_all(name, optim)
            % Plot the logged data across all iterations.
            
            % extract the data for all iterations (error metric and timing)
            conv_vec = 1:length(optim);
            for i=1:length(optim)
                err = optim{i}.sol_fom.err;
                t_solver = optim{i}.sol_fom.t_solver;
                
                if isfinite(err)
                    err_conv_ok_vec(i) = err;
                    err_conv_ko_vec(i) = NaN;
                else
                    err_conv_ok_vec(i) = NaN;
                    err_conv_ko_vec(i) = 0;
                end
                t_solver_conv_vec(i) = t_solver;
            end
                        
            % plot all iterations
            figure('name', [name ' / all'])
            
            % plot the error metric
            subplot(1,2,1)
            plot(conv_vec, 1e2.*err_conv_ok_vec, 'xg')
            hold('on')
            plot(conv_vec, 1e2.*err_conv_ko_vec, 'xr')
            grid('on')
            xlabel('i (#)')
            ylabel('err (%)')
            title(sprintf('Convergence / %s', name), 'interpreter', 'none')
            
            % plot the solver timing data
            subplot(1,2,2)
            plot(conv_vec, seconds(t_solver_conv_vec), 'xb')
            grid('on')
            xlabel('i (#)')
            ylabel('t (s)')
            title(sprintf('Time / %s', name), 'interpreter', 'none')
        end
    end
    
    %% private static api
    methods(Static, Access = private)
        function optim = get_log(x_unclamp, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter, data_optim)
            % Parse and assign the logging data for an iteration.
            
            % extract
            fct_err = data_optim.fct_err;
            fct_clamp = data_optim.fct_clamp;
            fct_param = data_optim.fct_param;
            n_var = data_optim.n_var;
            
            % transform unconstrained variables into bounded variables with sine transformation
            x_scale = fct_clamp(x_unclamp);
            pop_valid = isfinite(err);

            % select the best point
            [err, idx_best] = min(err);
            x_scale = x_scale(idx_best,:);

            % get the error metrics (handle empty/invalid points)
            if (isempty(x_scale)==true)||(isfinite(err)==false)
                x_scale = NaN(1, n_var);
                err_wgt_vec = NaN;
                err_vec = NaN;
                wgt_vec = NaN;
                err = NaN;
            else
                % get the error metrics
                [err, err_vec, wgt_vec] = fct_err(x_scale);
                                
                % get the weighted error vector
                err_wgt_vec = repelem(err_vec, round(wgt_vec));
            end
                        
            % extract the parameter structure from a raw data (transformation and normalization)
            [n_pts, param, bnd, is_bound] = fct_param(x_scale);
            assert(n_pts==1, 'invalid size: solution')
                        
            % parse the solver message
            msg = splitlines(strtrim(msg));
            
            % get the solver figures of merit
            sol_fom = SolverLog.get_sol_fom(is_valid, is_bound, err, n_iter, n_eval, pop_valid, t_solver, t_iter);
            
            % parse the error metrics
            err_fom = SolverLog.get_err_fom(err_vec, wgt_vec, err_wgt_vec);
            
            % assign
            optim.param = param;
            optim.bnd = bnd;
            optim.err_wgt_vec = err_wgt_vec;
            optim.err_vec = err_vec;
            optim.wgt_vec = wgt_vec;
            optim.err_fom = err_fom;
            optim.sol_fom = sol_fom;
            optim.msg = msg;
        end
        
        function sol_fom = get_sol_fom(is_valid, is_bound, err, n_iter, n_eval, pop_valid, t_solver, t_iter)
            % Get the solver figure of merit.
            
            % check if the point is valid
            is_valid = is_valid&&isfinite(err);

            % get the size of the solution (number of parameter combinations)
            n_pop_all = length(pop_valid);
            n_pop_fail = nnz(pop_valid==false);
            
            % assign
            sol_fom.err = err;
            sol_fom.is_bound = is_bound;
            sol_fom.is_valid = is_valid;
            sol_fom.n_iter = n_iter;
            sol_fom.n_eval = n_eval;
            sol_fom.n_pop_all = n_pop_all;
            sol_fom.n_pop_fail = n_pop_fail;
            sol_fom.t_solver = t_solver;
            sol_fom.t_iter = t_iter;
        end
        
        function err_fom = get_err_fom(err_vec, wgt_vec, err_wgt_vec)
            % Parse the error metrics.
            
            % get the error for different types of norms
            err_fom.norm_avg = SolverUtils.get_norm(err_vec, wgt_vec, 1);
            err_fom.norm_n_2 = SolverUtils.get_norm(err_vec, wgt_vec, 2);
            err_fom.norm_n_4 = SolverUtils.get_norm(err_vec, wgt_vec, 4);
            err_fom.norm_n_6 = SolverUtils.get_norm(err_vec, wgt_vec, 6);
            err_fom.norm_n_8 = SolverUtils.get_norm(err_vec, wgt_vec, 8);
            err_fom.norm_n_10 = SolverUtils.get_norm(err_vec, wgt_vec, 10);
            err_fom.norm_n_12 = SolverUtils.get_norm(err_vec, wgt_vec, 12);
            err_fom.norm_inf = SolverUtils.get_norm(err_vec, wgt_vec, Inf);
            
            % get error for different percentiles
            err_fom.percentile_50 = SolverUtils.get_percentile(err_wgt_vec, 0.50);
            err_fom.percentile_75 = SolverUtils.get_percentile(err_wgt_vec, 0.75);
            err_fom.percentile_90 = SolverUtils.get_percentile(err_wgt_vec, 0.90);
            err_fom.percentile_95 = SolverUtils.get_percentile(err_wgt_vec, 0.95);
            err_fom.percentile_99 = SolverUtils.get_percentile(err_wgt_vec, 0.99);
            
            % get the size of the error vector and the weighted error vector 
            err_fom.n_rep = length(err_wgt_vec);
            err_fom.n_all = length(err_vec);
        end
    end
end