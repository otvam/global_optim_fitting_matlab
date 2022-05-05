classdef SolverLog < handle
    %% properties
    properties (SetAccess = private, GetAccess = private)
        solver_type
        data_optim
        
        t_start
        t_last
        optim
    end
    
    %% public
    methods (Access = public)
        function self = SolverLog(solver_type, data_optim)
            % set data
            self.data_optim = data_optim;
            self.solver_type = solver_type;
            
            % init data
            self.t_start = datetime('now');
            self.t_last = datetime('now');
            self.optim = {};
            
        end
        
        function get_iter(self, x_unclamp, n_iter, n_eval, msg)
            % Get the initial values.
            
            t_now = datetime('now');
            t_solver = t_now-self.t_start;
            t_iter = t_now-self.t_last;
            self.t_last = t_now;
            
            is_valid = true;
            msg = sprintf('intermediate results\nstate: %s', msg);
            
            optim_tmp = SolverLog.get_log(x_unclamp, is_valid, n_iter, n_eval, msg, t_solver, t_iter, self.data_optim);
            self.optim{end+1} = optim_tmp;
            
            name = sprintf('iter / %s / n = %d / %d', self.solver_type, n_iter, n_eval);
            SolverLog.get_disp(name, self.optim{end})
        end
        
        function get_final(self, x_unclamp, n_iter, n_eval, msg, is_valid)
            % Get the initial values.
            
            t_now = datetime('now');
            
            t_solver = t_now-self.t_start;
            t_iter = t_solver./(n_iter+1);
            
            msg = sprintf('final results\n%s', msg);
            
            optim_tmp = SolverLog.get_log(x_unclamp, is_valid, n_iter, n_eval, msg, t_solver, t_iter, self.data_optim);
            self.optim{end+1} = optim_tmp;
            
            name = sprintf('final / %s', self.solver_type);
            SolverLog.get_disp(name, self.optim{end})
            SolverLog.get_plot_single(name, self.optim{end})
            SolverLog.get_plot_all(name, self.optim)
        end
        
        function optim = get_optim(self)
            % Get the log.
            
            optim = self.optim;
        end
    end
    
    %% public api
    methods(Static, Access = public)
        function get_disp(name, optim)
            % Display the optim structure.
            
            fprintf('    %s\n', name)
            
            fprintf('        msg\n')
            for i=1:length(optim.msg)
                fprintf('            %s\n', optim.msg{i})
            end
            
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
            
            fprintf('        param\n')
            field = fieldnames(optim.param);
            for i=1:length(field)
                value = optim.param.(field{i});
                
                if size(value, 2)==1
                    fprintf('            param.%s = %s;\n', field{i}, SolverUtils.get_disp_vec(value))
                else
                    fprintf('            param.%s = %s;\n', field{i}, SolverUtils.get_disp_nan(value))
                end
            end
            
            fprintf('        bnd\n')
            field = fieldnames(optim.bnd);
            for i=1:length(field)
                value = optim.bnd.(field{i});
                fprintf('            bnd.%s = %s;\n', field{i}, SolverUtils.get_disp_vec(value))
            end
        end
        
        function get_plot_single(name, optim)
            % Plot the optim structure.
            
            % extract last
            err_wgt_vec = optim.err_wgt_vec;
            err_vec = optim.err_vec;
                                    
            % plot last iteration
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
            % Plot the optim structure.
            
            % extract all
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
            
            subplot(1,2,1)
            semilogy(conv_vec, 1e2.*err_conv_ok_vec, 'xg')
            hold('on')
            semilogy(conv_vec, 1e2.*err_conv_ko_vec, 'xr')
            grid('on')
            xlabel('i (#)')
            ylabel('err (%)')
            title(sprintf('Convergence / %s', name), 'interpreter', 'none')
            
            subplot(1,2,2)
            plot(conv_vec, seconds(t_solver_conv_vec), 'xb')
            grid('on')
            xlabel('i (#)')
            ylabel('t (s)')
            title(sprintf('Time / %s', name), 'interpreter', 'none')
        end
    end
    
    %% private api
    methods(Static, Access = private)
        function optim = get_log(x_unclamp, is_valid, n_iter, n_eval, msg, t_solver, t_iter, data_optim)
            % Create the optim struct.
            
            % extract
            fct_opt = data_optim.fct_opt;
            fct_clamp = data_optim.fct_clamp;
            fct_param = data_optim.fct_param;
            error_norm = data_optim.error_norm;
            n_var = data_optim.n_var;
            
            % handle empty points
            if isempty(x_unclamp)==true
                x_unclamp = NaN(1, n_var);
            end
            
            % get error
            x_scale = fct_clamp(x_unclamp);
            [err_mat, wgt_mat] = fct_opt(x_scale);
            err = SolverUtils.get_norm(err_mat, wgt_mat, error_norm);
            pop_valid = isfinite(err);
            
            [err, idx] = min(err);
            x_scale = x_scale(idx,:); 
            err_vec = err_mat(idx,:); 
            wgt_vec = wgt_mat(idx,:); 
                        
            % handle invalid points
            if isfinite(err)==false
                x_scale = NaN(1, n_var);
                err_wgt_vec = NaN;
                err_vec = NaN;
                wgt_vec = NaN;
            else
                err_wgt_vec = repelem(err_vec, round(wgt_vec));
            end
                        
            % parse point
            [n_pts, param, bnd, is_bound] = fct_param(x_scale);
            assert(n_pts==1, 'invalid size: solution')
            
            % compute
            is_valid = is_valid&&isfinite(err);
            
            % parse message
            msg = splitlines(strtrim(msg));
            
            % compute
            sol_fom = SolverLog.get_sol_fom(is_valid, is_bound, err, n_iter, n_eval, pop_valid, t_solver, t_iter);
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
            % Get solution mectris.
            
            % compute
            is_valid = is_valid&&isfinite(err);
            
            % get size
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
            % Get error metrics.
            
            % norm error
            err_fom.norm_avg = SolverUtils.get_norm(err_vec, wgt_vec, 1);
            err_fom.norm_n_2 = SolverUtils.get_norm(err_vec, wgt_vec, 2);
            err_fom.norm_n_4 = SolverUtils.get_norm(err_vec, wgt_vec, 4);
            err_fom.norm_n_6 = SolverUtils.get_norm(err_vec, wgt_vec, 6);
            err_fom.norm_n_8 = SolverUtils.get_norm(err_vec, wgt_vec, 8);
            err_fom.norm_n_10 = SolverUtils.get_norm(err_vec, wgt_vec, 10);
            err_fom.norm_n_12 = SolverUtils.get_norm(err_vec, wgt_vec, 12);
            err_fom.norm_inf = SolverUtils.get_norm(err_vec, wgt_vec, Inf);
            
            % percentile error
            err_fom.percentile_50 = SolverUtils.get_percentile(err_wgt_vec, 0.50);
            err_fom.percentile_75 = SolverUtils.get_percentile(err_wgt_vec, 0.75);
            err_fom.percentile_90 = SolverUtils.get_percentile(err_wgt_vec, 0.90);
            err_fom.percentile_95 = SolverUtils.get_percentile(err_wgt_vec, 0.95);
            err_fom.percentile_99 = SolverUtils.get_percentile(err_wgt_vec, 0.99);
            
            % error size
            err_fom.n_rep = length(err_wgt_vec);
            err_fom.n_all = length(err_vec);
        end
    end
end