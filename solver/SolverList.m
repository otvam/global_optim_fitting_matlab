classdef SolverList < handle
    % Static class provided a common interface for different solvers.
    %
    %    Support fminunc / fminsearch / fmincon.
    %    Support surrogateopt / particleswarm / ga.
    %    Special init solver for finding reasonable initial values.
    %    Handle callback output function (display and logging).
    
    %    Thomas Guillod.
    %    2021-2022 - BSD License.
    
    %% public static api
    methods (Static, Access = public)
        function x = get_solver(fct_sol, fct_iter, fct_final, x0, lb, ub, options, solver_type)
            % Common interface for the different solvers.
            
            % call the specified solver
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
            
            % call the display and logging function with the final values
            fct_final(x, n_iter, n_eval, msg, is_valid);
        end
    end
    
    %% private static api
    methods (Static, Access = private)
        function [x, is_valid, n_iter, n_eval, msg] = get_init(fct_sol, fct_iter, x0, lb, ub, options)
            % Special dummy solver for finding several reasonable initial value.
            
            % solver need contrained variables
            n_init = size(x0, 1);
            n_var = size(x0, 2);
            is_bnd = all(isfinite(lb))&&all(isfinite(ub));
            assert(is_bnd, 'invalid initial point')
            assert(n_init>=1, 'invalid data')
            
            % extract
            n_batch = options.n_batch;
            n_tot = options.n_tot;
            n_eval_max = options.n_eval_max;
            val_lim = options.val_lim;
            
            % run solver
            n_iter = 0;
            n_eval = 0;
            x = [];
            err = [];
            while (n_eval<=n_eval_max)&&(size(x, 1)<n_tot)
                n_iter = n_iter+1;
                n_eval = n_eval+n_batch;
                
                [x, err] = SolverList.get_iter_init(x, err, n_iter, n_eval, fct_sol, fct_iter, n_var, lb, ub, n_batch, val_lim);
            end
            
            % check if the required number of solution is found
            if size(x, 1)==0
                x = x0;
                is_valid = false;
                msg = 'no solution found';
            elseif size(x, 1)<n_tot
                is_valid = false;
                msg = 'number of solutions is not sufficient';
            else
                [~, idx] = sort(err);
                x = x(idx,:);
                x = x(1:n_tot,:);

                is_valid = true;
                msg = 'number of solutions is sufficient';
            end
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_fminunc(fct_sol, fct_iter, x0, options)
            % Call the MATLAB fminunc solver.
            
            % solver need a single finite initial value
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            % call the solver
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, ~, exitflag, output] = fminunc(fct_sol, x0, options);
            
            % check convergence
            is_valid = any(exitflag==[1 2 3 5]);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_fminsearch(fct_sol, fct_iter, x0, options)
            % Call the MATLAB fminsearch solver.
            
            % solver need a single finite initial value
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            % call the solver
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, ~, exitflag, output] = fminsearch(fct_sol, x0, options);
            
            % check convergence
            is_valid = any(exitflag==1);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_fmincon(fct_sol, fct_iter, x0, lb, ub, options)
            % Call the MATLAB fmincon solver.
            
            % solver need a single finite initial value
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            % call the solver
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, ~, exitflag, output] = fmincon(fct_sol, x0, [], [], [], [], lb, ub, [], options);
            
            % check convergence
            is_valid = any(exitflag==[1 2 3 4 5]);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_surrogateopt(fct_sol, fct_iter, x0, lb, ub, options)
            % Call the MATLAB surrogateopt solver.
            
            % solver need contrained variables
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            is_bnd = all(isfinite(lb))&&all(isfinite(ub));
            assert(is_bnd, 'invalid initial point')
            assert(n_init>=1, 'invalid data')
            
            % call the solver (assign initial values if possible)
            options = optimoptions(options, 'Display', 'off');
            options = optimoptions(options, 'PlotFcn', []);
            options = optimoptions(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            if is_init==true
                options = optimoptions(options, 'InitialPoints', x0);
            end
            [x, ~, exitflag, output] = surrogateopt(fct_sol, lb, ub, [], [], [], [], [], options);
            
            % check convergence
            is_valid = any(exitflag==[1 3 10]);
            n_iter = output.funccount;
            n_eval = output.funccount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_particleswarm(fct_sol, fct_iter, x0, lb, ub, options)
            % Call the MATLAB particleswarm solver.
            
            % solver has no particular constraints
            n_var = size(x0, 2);
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(n_init>=1, 'invalid data')
            
            % call the solver (assign initial values if possible)
            options = optimoptions(options, 'Display', 'off');
            options = optimoptions(options, 'PlotFcn', []);
            options = optimoptions(options, 'OutputFcn', @(optim, state) SolverList.get_outfun_particleswarm(optim, state, fct_iter));
            if is_init==true
                options = optimoptions(options, 'InitialSwarmMatrix', x0);
            end
            [x, ~, exitflag, output] = particleswarm(fct_sol, n_var, lb, ub, options);
            
            % check convergence
            is_valid = any(exitflag==1);
            n_iter = output.iterations;
            n_eval = output.funccount;
            msg = output.message;
        end
        
        function [x, is_valid, n_iter, n_eval, msg] = get_ga(fct_sol, fct_iter, x0, lb, ub, options)
            % Call the MATLAB ga solver.
            
            % solver has no particular constraints
            n_var = size(x0, 2);
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(n_init>=1, 'invalid data')
            
            % call the solver (assign initial values if possible)
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
            % Output function for standard MATLAB solvers.
            
            % extract
            n_iter = optim.iteration;
            n_eval = optim.funccount;
            
            % call the display and logging function
            fct_iter(x, n_iter, n_eval, msg);
            
            % do not stop the solver
            stop = false;
        end
        
        function stop = get_outfun_particleswarm(optim, msg, fct_iter)
            % Output function for particleswarm.
            
            % extract
            n_iter = optim.iteration;
            n_eval = optim.funccount;
            x = optim.bestx;
            
            % call the display and logging function
            fct_iter(x, n_iter, n_eval, msg);
            
            % do not stop the solver
            stop = false;
        end
        
        function [optim, options, optchanged] = get_outfun_ga(options, optim, msg, fct_iter)
            % Output function for ga.
            
            % extract
            n_iter = optim.Generation;
            n_eval = optim.FunEval;
            x = optim.Population;
            
            % call the display and logging function
            fct_iter(x, n_iter, n_eval, msg);
            
            % do not stop the solver
            optchanged = false;
        end
        
        function [x, err] = get_iter_init(x, err, n_iter, n_eval, fct_sol, fct_iter, n_var, lb, ub, n_batch, val_lim)
            % Make an iteration for the init solver.
            
            % select random points between the bounds
            x_tmp = lb+(ub-lb).*rand(n_batch, n_var);
            
            % evaluate the error function
            err_tmp = fct_sol(x_tmp);
            
            % select only points that are valid a below the threshold
            idx = isfinite(err_tmp)&(err_tmp<val_lim);
            x = [x ; x_tmp(idx,:)];
            err = [err ; err_tmp(idx,:)];
            
            % call the display and logging function
            fct_iter(x, n_iter, n_eval, 'iter');
        end
    end
end