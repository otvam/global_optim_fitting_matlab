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
        log_iter % log (or not) the solver iterations
        log_final % log (or not) the solver final results
        fct_unscale % function the extract the parameter structure from a raw matrix
        fct_sol % function returning the error metrics
        format % structure with formatting instructions (name and unit)
        
        t_start % timestamp set at the initialization
        t_last % timestamp of the last iteration
        optim % cell containing the logged data
    end
    
    %% public
    methods (Access = public)
        function self = SolverLog(solver_type, log_iter, log_final, fct_unscale, fct_sol, format)
            % Constructor.
            
            % set data
            self.solver_type = solver_type;
            self.log_iter = log_iter;
            self.log_final = log_final;
            self.fct_unscale = fct_unscale;
            self.fct_sol = fct_sol;
            self.format = format;
            
            % init the timing and logging data
            self.t_start = datetime('now');
            self.t_last = datetime('now');
            self.optim = {};
        end
        
        function get_iter(self, x_scale, err, n_iter, n_eval, msg, is_valid)
            % Log and display the solver progress after an iteration.
            
            % get if logging is required
            if self.log_iter==false
                return
            end
            
            % get the total solver timing
             [t_solver, t_iter] = get_time(self, n_iter, false);
                         
            % get message
            msg = sprintf('intermediate results\nstate: %s', msg);
            
            % get and add the logging data
            optim_tmp = self.get_log(x_scale, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter);
            self.optim{end+1} = optim_tmp;
            
            % display the data
            name = sprintf('iter / %s / n = %d / %d', self.solver_type, n_iter, n_eval);
            SolverLog.get_disp(name, self.optim{end}, self.format)
        end
        
        function get_final(self, x_scale, err, n_iter, n_eval, msg, is_valid)
            % Log, display, and plot the solver results after the final iteration.
            
            % get if logging is required
            if self.log_final==false
                return
            end
            
            % get the total solver timing
             [t_solver, t_iter] = get_time(self, n_iter, true);
            
            % get message
            msg = sprintf('final results\n%s', msg);
            
            % get and add the logging data
            optim_tmp = self.get_log(x_scale, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter);
            self.optim{end+1} = optim_tmp;
            
            % display and plot the data
            name = sprintf('final / %s', self.solver_type);
            SolverLog.get_disp(name, self.optim{end}, self.format)
            SolverLog.get_plot_single(name, self.optim{end}, self.format)
            SolverLog.get_plot_all(name, self.optim, self.format)
        end
        
        function optim = get_optim(self)
            % Get the logged data.
            
            optim = self.optim;
        end
    end
    
    %% public static api
    methods(Static, Access = public)
        function get_disp(name, optim, format)
            % Display the logged data for a specific iteration.
            
            % name of the iteration
            fprintf('    %s\n', name)
            err = format.err;
            param = format.param;
            
            % display the solver figures of merit
            fprintf('        fom\n')
            fprintf('            msg\n')
            for i=1:length(optim.sol_fom.msg)
                fprintf('                %s\n', optim.sol_fom.msg{i})
            end
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
            fprintf('                n_pop_all = %d\n', optim.sol_fom.n_pop_all)
            fprintf('                n_pop_fail = %d\n', optim.sol_fom.n_pop_fail)
            fprintf('                err_best = %s\n', SolverLog.get_format_scalar(optim.sol_fom.err_best, err))
            
            if optim.has_solution==true
                % display the error metrics
                fprintf('        err_fom\n')
                fprintf('            size\n')
                fprintf('                n_set = %d\n', optim.err_fom.n_set)
                fprintf('                err_best = %s\n', SolverLog.get_format_scalar(optim.sol_fom.err_best, err))
                fprintf('            weight\n')
                fprintf('                wgt_sum = %.3f\n', optim.err_fom.wgt_sum)
                fprintf('                wgt_avg = %.3f\n', optim.err_fom.wgt_avg)
                if optim.err_fom.n_set>1
                    fprintf('            avg\n')
                    fprintf('                avg = %s\n', SolverLog.get_format_scalar(optim.err_fom.avg, err))
                    fprintf('                min = %s\n', SolverLog.get_format_scalar(optim.err_fom.min, err))
                    fprintf('                max = %s\n', SolverLog.get_format_scalar(optim.err_fom.max, err))
                    fprintf('            norm\n')
                    for i=1:length(optim.err_fom.norm_val)
                        str_val = sprintf('%d', optim.err_fom.norm_val(i));
                        str_err = SolverLog.get_format_scalar(optim.err_fom.norm_err(i), err);
                        fprintf('                norm / %s = %s\n', str_val, str_err)
                    end
                    fprintf('            percentile\n')
                    for i=1:length(optim.err_fom.percentile_val)
                        str_val = sprintf('%.1f %%', 1e2.*optim.err_fom.percentile_val(i));
                        str_err = SolverLog.get_format_scalar(optim.err_fom.percentile_err(i), err);
                        fprintf('                percentile / %s = %s\n', str_val, str_err)
                    end
                end
                
                % display the current best parameter combinations
                fprintf('        param\n')
                field = fieldnames(optim.param);
                for i=1:length(field)
                    value = optim.param.(field{i});
                    
                    if isfield(param, field{i})
                        str = SolverLog.get_format_vec(value, param.(field{i}));
                        fprintf('            param.%s = %s\n', field{i}, str)
                    else
                        fprintf('            param.%s = hidden\n', field{i})
                    end
                end
                
                % display if the parameters are close to the bounds
                fprintf('        bnd\n')
                field = fieldnames(optim.bnd);
                for i=1:length(field)
                    value = optim.bnd.(field{i});
                    str = SolverLog.get_format_bnd(value);
                    fprintf('            bnd.%s = %s\n', field{i}, str)
                end
            end
        end
        
        function get_plot_single(name, optim, format)
            % Plot the logged data for a specific iteration.
            
            if optim.has_solution==true
                % extract format
                scale = format.err.scale;
                unit = format.err.unit;
                
                % extract the error vector and the weighted error vector
                err_vec = optim.err_fom.err_vec;
                wgt_vec = optim.err_fom.wgt_vec;
                n_set = optim.err_fom.n_set;
                
                % plot the error distribution (if it exists)
                if n_set>1
                    figure('name', [name ' / single'])
                    
                    % histogram
                    subplot(1,2,1)
                    histogram(scale.*err_vec, 'Normalization', 'pdf')
                    grid('on')
                    xlabel(['err (' unit ')'], 'interpreter', 'none')
                    ylabel('p (#)', 'interpreter', 'none')
                    title(sprintf('Histogram / %s', name), 'interpreter', 'none')
                    
                    % scatter plot with error and weights
                    subplot(1,2,2)
                    pts_vec = 1:n_set;
                    [err_vec, idx] = sort(err_vec);
                    wgt_vec = wgt_vec(idx);
                    scatter(pts_vec, scale.*err_vec, 50, wgt_vec, 'filled')
                    colorbar();
                    xlim([min(pts_vec) max(pts_vec)])
                    grid('on')
                    xlabel('idx (#)', 'interpreter', 'none')
                    ylabel(['err (' unit ')'], 'interpreter', 'none')
                    title(sprintf('Weight / %s', name), 'interpreter', 'none')
                end
            end
        end
        
        function get_plot_all(name, optim, format)
            % Plot the logged data across all iterations.
            
            % extract format
            scale = format.err.scale;
            unit = format.err.unit;
            
            % extract the data for all iterations (error metric and timing)
            conv_vec = 1:length(optim);
            for i=1:length(optim)
                err_best = optim{i}.sol_fom.err_best;
                t_solver = optim{i}.sol_fom.t_solver;
                
                err_conv_vec(i) = err_best;
                t_solver_conv_vec(i) = t_solver;
            end
            
            % plot all iterations
            figure('name', [name ' / all'])
            
            % plot the error metric
            subplot(1,2,1)
            plot(conv_vec, scale.*err_conv_vec, 'og')
            xlim([min(conv_vec) max(conv_vec)])
            grid('on')
            xlabel('iter (#)', 'interpreter', 'none')
            ylabel(['err (' unit ')'], 'interpreter', 'none')
            title(sprintf('Convergence / %s', name), 'interpreter', 'none')
            
            % plot the solver timing data
            subplot(1,2,2)
            plot(conv_vec, seconds(t_solver_conv_vec), 'ob')
            xlim([min(conv_vec) max(conv_vec)])
            grid('on')
            xlabel('iter (#)', 'interpreter', 'none')
            ylabel('t (s)', 'interpreter', 'none')
            title(sprintf('Time / %s', name), 'interpreter', 'none')
        end
    end
    
    %% private api
    methods( Access = private)
        function [t_solver, t_iter] = get_time(self, n_iter, is_final)
            % Compute and update the solver timing (total time and iteration time).
            
            % get current time
            t_now = datetime('now');
            
            % get total solver time
            t_solver = t_now-self.t_start;
            
            % the time of the last iteration or average time per iteration
            if is_final==true
                t_iter = t_solver./(n_iter+1);
            else
                t_iter = t_now-self.t_last;
            end
            
            % update the iteration timestamp
            self.t_last = t_now;
        end
        
        function optim = get_log(self, x_scale, err, is_valid, n_iter, n_eval, msg, t_solver, t_iter)
            % Parse and assign the logging data for a specific iteration.
            
            % check if there is a valid solution
            has_solution = (isempty(x_scale)==false)&&(isempty(err)==false)&&any(isfinite(err));
            
            % parse the best solution (if any)
            if has_solution==true
                % best the best parameter combination
                [~, idx_best] = min(err);
                x_scale = x_scale(:,idx_best);
                
                % extract the parameter structure from a raw data
                [n_pts, param, bnd, is_bound] = self.fct_unscale(x_scale);
                assert(n_pts==1, 'invalid size: solution')
                
                % get the error metrics
                [err_best, err_vec, wgt_vec] = self.fct_sol(x_scale);
                
                % process the error metrics
                err_fom = SolverLog.get_err_fom(err_best, err_vec, wgt_vec);
            else
                param = [];
                bnd = [];
                err_fom = [];
                err_best = NaN;
                is_valid = false;
                is_bound = false;
            end
            
            % get the solver figures of merit
            sol_fom = SolverLog.get_sol_fom(err, err_best, is_valid, is_bound, n_iter, n_eval, msg, t_solver, t_iter);
            
            % assign base data
            optim.has_solution = has_solution;
            optim.sol_fom = sol_fom;
            optim.param = param;
            optim.bnd = bnd;
            optim.err_fom = err_fom;
        end
    end
    
    %% private static api
    methods(Static, Access = private)
        function txt = get_format_scalar(val, format)
            % Parse a scalar to a string with scaling and units.
            
            % extract
            spec = format.spec;
            scale = format.scale;
            unit = format.unit;
            
            % parse each elements
            txt = sprintf(spec, scale.*val);
            
            % add unit
            txt = [txt ' ' unit];
        end
        
        function txt = get_format_vec(vec, format)
            % Parse a vector to a string with scaling and units.
            
            % extract
            spec = format.spec;
            scale = format.scale;
            unit = format.unit;
            
            % parse each elements
            for i=1:length(vec)
                txt{i} = sprintf(spec, scale.*vec(i));
            end
            
            % assemble the string
            txt = sprintf('[%s]', strjoin(txt, ' ; '));
            
            % add unit
            txt = [txt ' ' unit];
        end
        
        function txt = get_format_bnd(vec)
            % Parse a boolean vector to string.
            
            % parse each elements
            for i=1:length(vec)
                txt{i} = mat2str(vec(i));
            end
            
            % assemble the string
            txt = sprintf('[%s]', strjoin(txt, ' ; '));
        end
        
        function sol_fom = get_sol_fom(err, err_best, is_valid, is_bound, n_iter, n_eval, msg, t_solver, t_iter)
            % Process the solver data and assign the results to a struct.

            % check population
            pop_valid = isfinite(err);
            n_pop_all = length(pop_valid);
            n_pop_fail = nnz(pop_valid==false);
            
            % parse the solver message
            msg = splitlines(strtrim(msg));
            
            % get the solver figures of merit
            sol_fom.msg = msg;
            sol_fom.err_best = err_best;
            sol_fom.is_bound = is_bound;
            sol_fom.is_valid = is_valid;
            sol_fom.n_iter = n_iter;
            sol_fom.n_eval = n_eval;
            sol_fom.n_pop_all = n_pop_all;
            sol_fom.n_pop_fail = n_pop_fail;
            sol_fom.t_solver = t_solver;
            sol_fom.t_iter = t_iter;
        end
        
        function err_fom = get_err_fom(err_best, err_vec, wgt_vec)
            % Process the error metrics and assign the results to a struct.
            
            % get the error for different types of norms
            err_fom.n_set = (length(err_vec)+length(wgt_vec))./2;
            err_fom.wgt_sum = sum(wgt_vec);
            err_fom.wgt_avg = sum(wgt_vec)./length(wgt_vec);
            
            if err_fom.n_set>1
                % get error with simple metrics
                err_fom.avg = SolverUtils.get_error(err_vec, wgt_vec, 'avg');
                err_fom.min = SolverUtils.get_error(err_vec, wgt_vec, 'min');
                err_fom.max = SolverUtils.get_error(err_vec, wgt_vec, 'max');
                
                % get error for different norm
                norm = [1 2 4 6 8 10 12 Inf];
                for i=1:length(norm)
                    err_fom.norm_err(i) = SolverUtils.get_norm(err_vec, wgt_vec, norm(i));
                    err_fom.norm_val(i) = norm(i);
                end
                
                % get error for different percentiles
                percentile = [0.01 0.05 0.1 0.5 0.9 0.95 0.99];
                for i=1:length(percentile)
                    err_fom.percentile_err(i) = SolverUtils.get_percentile(err_vec, wgt_vec, percentile(i));
                    err_fom.percentile_val(i) = percentile(i);
                end
            end
            
            % assign raw data
            err_fom.err_best = err_best;
            err_fom.err_vec = err_vec;
            err_fom.wgt_vec = wgt_vec;
        end
    end
end