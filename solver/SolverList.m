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
                    [x, err, is_valid, n_iter, n_eval, msg] = SolverList.get_init(fct_sol, fct_iter, x0, lb, ub, options);
                case 'fminunc'
                    [x, err, is_valid, n_iter, n_eval, msg] = SolverList.get_fminunc(fct_sol, fct_iter, x0, options);
                case 'fminsearch'
                    [x, err, is_valid, n_iter, n_eval, msg] = SolverList.get_fminsearch(fct_sol, fct_iter, x0, options);
                case 'fmincon'
                    [x, err, is_valid, n_iter, n_eval, msg] = SolverList.get_fmincon(fct_sol, fct_iter, x0, lb, ub, options);
                case 'surrogateopt'
                    [x, err, is_valid, n_iter, n_eval, msg] = SolverList.get_surrogateopt(fct_sol, fct_iter, x0, lb, ub, options);
                case 'particleswarm'
                    [x, err, is_valid, n_iter, n_eval, msg] = SolverList.get_particleswarm(fct_sol, fct_iter, x0, lb, ub, options);
                case 'ga'
                    [x, err, is_valid, n_iter, n_eval, msg] = SolverList.get_ga(fct_sol, fct_iter, x0, lb, ub, options);
                otherwise
                    error('invalid data')
            end
            
            % call the display and logging function with the final values
            fct_final(x, err, n_iter, n_eval, msg, is_valid);
        end
    end
    
    %% private static api
    methods (Static, Access = private)
        function [x, err, is_valid, n_iter, n_eval, msg] = get_init(fct_sol, fct_iter, x0, lb, ub, options)
            % Special dummy solver for finding several reasonable initial value.
            
            % solver need contrained variables
            n_init = size(x0, 1);
            n_var = size(x0, 2);
            is_init = all(isfinite(x0(:)));
            is_bnd = all(isfinite(lb))&&all(isfinite(ub));
            assert(is_bnd, 'invalid initial point')
            assert(n_init>=1, 'invalid data')
            
            % extract
            n_batch = options.n_batch;
            n_tot = options.n_tot;
            n_iter_max = options.n_iter_max;
            err_lim = options.err_lim;
            
            % run solver
            n_iter = 0;
            n_eval = 0;
            x = [];
            err = [];
            is_timeout = false;
            while (is_timeout==false)&&(n_iter<=n_iter_max)&&(size(x, 1)<n_tot)
                % select the point
                if (n_iter==0)&&(is_init==true)
                    x_tmp = x0;
                else
                    x_tmp = lb+(ub-lb).*rand(n_batch, n_var);
                end
                
                % eval the error function
                [is_timeout, x, err, n_eval] = SolverList.get_iter_init(x, err, x_tmp, n_iter, n_eval, fct_sol, fct_iter, err_lim);
                
                % update iteration
                n_iter = n_iter+1;
            end
            
            % check if the required number of solution is found
            if size(x, 1)==0
                x = x0;
                err = fct_sol(x0);
                is_valid = false;
                msg = 'no solution found';
            elseif size(x, 1)<n_tot
                is_valid = false;
                msg = 'number of solutions is not sufficient';
            else
                [err, idx] = sort(err);
                x = x(idx,:);
                x = x(1:n_tot,:);
                err = err(1:n_tot,:);
                
                is_valid = true;
                msg = 'number of solutions is sufficient';
            end
        end
        
        function [x, err, is_valid, n_iter, n_eval, msg] = get_fminunc(fct_sol, fct_iter, x0, options)
            % Call the MATLAB fminunc solver.
            
            % solver need a single finite initial value
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            % call the solver
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, err, exitflag, output] = fminunc(fct_sol, x0, options);
            
            % check convergence
            is_valid = any(exitflag==[1 2 3 5]);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, err, is_valid, n_iter, n_eval, msg] = get_fminsearch(fct_sol, fct_iter, x0, options)
            % Call the MATLAB fminsearch solver.
            
            % solver need a single finite initial value
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            % call the solver
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, err, exitflag, output] = fminsearch(fct_sol, x0, options);
            
            % check convergence
            is_valid = any(exitflag==1);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, err, is_valid, n_iter, n_eval, msg] = get_fmincon(fct_sol, fct_iter, x0, lb, ub, options)
            % Call the MATLAB fmincon solver.
            
            % solver need a single finite initial value
            n_init = size(x0, 1);
            is_init = all(isfinite(x0(:)));
            assert(is_init, 'invalid initial point')
            assert(n_init==1, 'invalid data')
            
            % call the solver
            options = optimset(options, 'Display', 'off');
            options = optimset(options, 'OutputFcn', @(x, optim, state) SolverList.get_outfun_grad(x, optim, state, fct_iter));
            [x, err, exitflag, output] = fmincon(fct_sol, x0, [], [], [], [], lb, ub, [], options);
            
            % check convergence
            is_valid = any(exitflag==[1 2 3 4 5]);
            n_iter = output.iterations;
            n_eval = output.funcCount;
            msg = output.message;
        end
        
        function [x, err, is_valid, n_iter, n_eval, msg] = get_surrogateopt(fct_sol, fct_iter, x0, lb, ub, options)
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
            [x, err, exitflag, output] = surrogateopt(fct_sol, lb, ub, [], [], [], [], [], options);
            
            % check convergence
            is_valid = any(exitflag==[1 3 10]);
            n_iter = output.funccount;
            n_eval = output.funccount;
            msg = output.message;
        end
        
        function [x, err, is_valid, n_iter, n_eval, msg] = get_particleswarm(fct_sol, fct_iter, x0, lb, ub, options)
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
            [x, err, exitflag, output] = particleswarm(fct_sol, n_var, lb, ub, options);
            
            % check convergence
            is_valid = any(exitflag==1);
            n_iter = output.iterations;
            n_eval = output.funccount;
            msg = output.message;
        end
        
        function [x, err, is_valid, n_iter, n_eval, msg] = get_ga(fct_sol, fct_iter, x0, lb, ub, options)
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
            [x, err, exitflag, output] = ga(fct_sol, n_var, [], [], [], [], lb, ub, [], [], options);
            
            is_valid = any(exitflag==[1 3 4 5]);
            n_iter = output.generations;
            n_eval = output.funccount;
            msg = output.message;
        end
    end
    
    methods (Static, Access = private)
        function is_timeout = get_outfun_grad(x, optim, msg, fct_iter)
            % Output function for standard MATLAB solvers.
            
            % extract
            n_iter = optim.iteration;
            n_eval = optim.funccount;
            err = optim.fval;
            
            % call the display and logging function
            is_valid = isempty(x)==false;
            is_timeout = fct_iter(x, err, n_iter, n_eval, msg, is_valid);
        end
        
        function is_timeout = get_outfun_particleswarm(optim, msg, fct_iter)
            % Output function for particleswarm.
            
            % extract
            n_iter = optim.iteration;
            n_eval = optim.funccount;
            x = [optim.bestx ; optim.swarm];
            err = [optim.bestfval ; optim.swarmfvals];
            
            % call the display and logging function
            is_valid = isempty(x)==false;
            is_timeout = fct_iter(x, err, n_iter, n_eval, msg, is_valid);
        end
        
        function [optim, options, optchanged] = get_outfun_ga(options, optim, msg, fct_iter)
            % Output function for ga.
            
            % extract
            n_iter = optim.Generation;
            n_eval = optim.FunEval;
            x = optim.Population;
            err = optim.Score;
            
            % call the display and logging function
            is_valid = isempty(x)==false;
            is_timeout = fct_iter(x, err, n_iter, n_eval, msg, is_valid);
            
            % stop the solver if required
            if is_timeout==true
                optim.StopFlag = 'timeout';
            end
            
            % do not stop the solver options
            optchanged = false;
        end
        
        function [is_timeout, x, err, n_eval] = get_iter_init(x, err, x_tmp, n_iter, n_eval, fct_sol, fct_iter, err_lim)
            % Make an iteration for the init solver.
            
            % evaluate the error function
            err_tmp = fct_sol(x_tmp);
            
            % count evaluation
            n_eval = n_eval+size(x_tmp, 1);
            
            % select only points that are valid a below the threshold
            idx = isfinite(err_tmp)&(err_tmp<err_lim);
            x = [x ; x_tmp(idx,:)];
            err = [err ; err_tmp(idx,:)];
            
            % call the display and logging function
            msg = 'iter';
            is_valid = isempty(x)==false;
            is_timeout = fct_iter(x, err, n_iter, n_eval, msg, is_valid);
        end
    end
end