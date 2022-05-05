classdef SolverList < handle
    %% public api
    methods (Static, Access = public)
        function x = get_solver(fct_sol, fct_iter, fct_final, x0, lb, ub, options, solver_type)
            switch solver_type
                case 'init'
                    [x, is_valid, n_iter, n_eval, msg] = SolverList.get_init(fct_sol, fct_iter, x0, lb, ub, options);
                case 'fminunc'
                    [x, is_valid, n_iter, n_eval, msg] = SolverList.get_fminunc(fct_sol, fct_iter, x0, options);
                case 'fminsearch'
                    [x, is_valid, n_iter, n_eval, msg] = SolverList.get_fminsearch(fct_sol, fct_iter, x0, options);
                case 'fmincon'
                    [x, is_valid, n_iter, n_eval, msg] = SolverList.get_fmincon(fct_sol, fct_iter, x0, lb, ub, options);
                case 'surrogateopt'
                    [x, is_valid, n_iter, n_eval, msg] = SolverList.get_surrogateopt(fct_sol, fct_iter, x0, lb, ub, options);
                case 'particleswarm'
                    [x, is_valid, n_iter, n_eval, msg] = SolverList.get_particleswarm(fct_sol, fct_iter, x0, lb, ub, options);
                case 'ga'
                    [x, is_valid, n_iter, n_eval, msg] = SolverList.get_ga(fct_sol, fct_iter, x0, lb, ub, options);
                otherwise
                    error('invalid data')
            end
            
            
            fct_final(x, n_iter, n_eval, msg, is_valid);
        end
    end
    
    %% private api
    methods (Static, Access = private)
        function [x, is_valid, n_iter, n_eval, msg] = get_init(fct_sol, fct_iter, x0, lb, ub, options)
            % check initial
            n_init = size(x0, 1);
            n_var = size(x0, 2);
            is_bnd = all(isfinite(lb))&&all(isfinite(ub));
            
            % extract
            n_batch = options.n_batch;
            n_tot = options.n_tot;
            n_eval_max = options.n_eval_max;
            val_lim = options.val_lim;
            
            % run solver
            assert(is_bnd, 'invalid initial point')
            assert(n_init>=1, 'invalid data')

            n_iter = 0;
            n_eval = 0;
            x = [];
            while (n_eval<=n_eval_max)&&(size(x, 1)<n_tot)
                n_iter = n_iter+1;
                n_eval = n_eval+n_batch;

                x = SolverList.get_iter_init(x, n_iter, n_eval, fct_sol, fct_iter, n_var, lb, ub, n_batch, val_lim);
            end
            
            if size(x, 1)==0
                x = x0;
                is_valid = false;
                msg = 'no solution found';
            elseif size(x, 1)<n_tot
                is_valid = false;
                msg = 'number of solutions is not sufficient';
            else                
                is_valid = true;
                msg = 'number of solutions is sufficient';
            end
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_fminunc(fct_sol, fct_iter, x0, options)
            % check initial
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            
            % run solver
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, ~, exitflag, output] = fminunc(fct_sol, x0, options);
            
            is_valid = any(exitflag==[1 2 3 5]);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_fminsearch(fct_sol, fct_iter, x0, options)
            % check initial
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            
            % run solver
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, ~, exitflag, output] = fminsearch(fct_sol, x0, options);
            
            is_valid = any(exitflag==1);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_fmincon(fct_sol, fct_iter, x0, lb, ub, options)
            % check initial
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            
            % run solver
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, ~, exitflag, output] = fmincon(fct_sol, x0, [], [], [], [], lb, ub, [], options);
            
            is_valid = any(exitflag==[1 2 3 4 5]);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_surrogateopt(fct_sol, fct_iter, x0, lb, ub, options)
            % check initial
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            
            % run solver
            assert(n_init>=1, 'invalid data')
            
            options = optimoptions(options, 'Display', 'off');
            options = optimoptions(options, 'PlotFcn', []);
            options = optimoptions(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            if is_init==true
                options = optimoptions(options, 'InitialPoints', x0);
            end
            [x, ~, exitflag, output] = surrogateopt(fct_sol, lb, ub, [], [], [], [], [], options);
            
            is_valid = any(exitflag==[1 3 10]);
            n_iter = output.funccount;
            n_eval = output.funccount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_particleswarm(fct_sol, fct_iter, x0, lb, ub, options)
            % check initial
            n_var = size(x0, 2);
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            
            % run solver
            assert(n_init>=1, 'invalid data')
            
            options = optimoptions(options, 'Display', 'off');
            options = optimoptions(options, 'PlotFcn', []);
            options = optimoptions(options, 'OutputFcn', @(optim, state) SolverList.get_outfun_particleswarm(optim, state, fct_iter));
            if is_init==true
                options = optimoptions(options, 'InitialSwarmMatrix', x0);
            end
            [x, ~, exitflag, output] = particleswarm(fct_sol, n_var, lb, ub, options);
            
            is_valid = any(exitflag==1);
            n_iter = output.iterations;
            n_eval = output.funccount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_ga(fct_sol, fct_iter, x0, lb, ub, options)
            % check initial
            n_var = size(x0, 2);
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            
            % run solver
            assert(n_init>=1, 'invalid data')
            
            options = optimoptions(options, 'Display', 'off');
            options = optimoptions(options, 'PlotFcn', []);
            options = optimoptions(options, 'OutputFcn', @(options, optim, state) SolverList.get_outfun_ga(options, optim, state, fct_iter));
            if is_init==true
                options = optimoptions(options, 'InitialPopulationMatrix', x0);
            end
            [x, ~, exitflag, output] = ga(fct_sol, n_var, [], [], [], [], lb, ub, [], [], options);
            
            is_valid = any(exitflag==[1 3 4 5]);
            n_iter = output.generations;
            n_eval = output.funccount;
            msg = output.message;
        end
    end
    
    methods (Static, Access = private)
        function stop = get_outfun_grad(x, optim, msg, fct_iter)
            % Output function for standard optimizers.
            
            n_iter = optim.iteration;
            n_eval = optim.funccount;
            
            fct_iter(x, n_iter, n_eval, msg);
            
            stop = false;
        end
        
        function stop = get_outfun_particleswarm(optim, msg, fct_iter)
            % Output function for particleswarm.
            
            n_iter = optim.iteration;
            n_eval = optim.funccount;
            x = optim.bestx;
            
            fct_iter(x, n_iter, n_eval, msg);
            
            stop = false;
        end
        
        function [optim, options, optchanged] = get_outfun_ga(options, optim, msg, fct_iter)
            % Output function for ga.
            
            n_iter = optim.Generation;
            n_eval = optim.FunEval;
            x = optim.Population;
                        
            fct_iter(x, n_iter, n_eval, msg);
            
            optchanged = false;
        end
        
        function x = get_iter_init(x, n_iter, n_eval, fct_sol, fct_iter, n_var, lb, ub, n_batch, val_lim)
            % Make an iteration.

            x_tmp = lb+(ub-lb).*rand(n_batch, n_var);
            err_tmp = fct_sol(x_tmp);

            idx = isfinite(err_tmp)&(err_tmp<val_lim);
            x = [x ; x_tmp(idx,:)];
                        
            msg = 'iter';
            fct_iter(x, n_iter, n_eval, msg);
        end
    end
end