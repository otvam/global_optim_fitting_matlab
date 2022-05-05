classdef SolverVar < handle
    %% properties
    properties (SetAccess = private, GetAccess = private)
        var_opt
        var_fix
        tol_bound
        
        x0_scale
        lb_scale
        ub_scale
    end
    
    %% public
    methods (Access = public)
        function self = SolverVar(var_opt, var_fix, tol_bound)
            % set data
            self.var_opt = var_opt;
            self.var_fix = var_fix;
            self.tol_bound = tol_bound;
            
            % init data
            for i=1:length(self.var_opt)
                [self.x0_scale(:,i), self.lb_scale(i), self.ub_scale(i)] = SolverVar.get_init(self.var_opt{i});
            end
        end
        
        function x0_scale = get_x0_scale(self)
            % Get the initial values.
            
            x0_scale = self.x0_scale;
        end
        
        function [n_pts, param, bnd, is_bound] = get_param(self, x_scale)
            % Get the parameters and bounds.
                
            for i=1:length(self.var_opt)
                [name, idx, x_tmp, is_bound_tmp] = SolverVar.get_param_opt(x_scale(:,i), self.var_opt{i}, self.lb_scale(i), self.ub_scale(i), self.tol_bound);
                                
                param.(name)(idx,:) = x_tmp;
                bnd.(name)(idx,:) = is_bound_tmp;
                is_bound(i,:) = is_bound_tmp;
            end
            
            n_rep = size(x_scale, 1);
            for i=1:length(self.var_fix)
                [name, idx, x_tmp, is_bound_tmp] = SolverVar.get_param_fix(self.var_fix{i}, n_rep);
                
                param.(name)(idx,:) = x_tmp;
                bnd.(name)(idx,:) = is_bound_tmp;
            end
            
            % order (otherwise random)
            n_pts = size(x_scale, 1);
            param = orderfields(param);
            bnd = orderfields(bnd);
            is_bound = all(is_bound, 1);
        end
        
        function [x_unclamp, lb_unclamp, ub_unclamp] = get_unclamp(self, x_scale, clamp_bnd)
            % Unclamp the variables.
            
            for i=1:size(x_scale, 2)
                [x_unclamp(:,i), lb_unclamp(i), ub_unclamp(i)] = SolverUtils.get_var_unclamp(x_scale(:,i), clamp_bnd, self.lb_scale(i), self.ub_scale(i));
            end
        end
        
        function x_scale = get_clamp(self, x_unclamp, clamp_bnd)
            % Clamp the variables.
                        
            for i=1:size(x_unclamp, 2)
                x_scale(:,i) = SolverUtils.get_var_clamp(x_unclamp(:,i), clamp_bnd, self.lb_scale(i), self.ub_scale(i));
            end
        end
    end
    
    %% private api
    methods(Static, Access = private)
        function [x0_scale, lb_scale, ub_scale] = get_init(var)
            % Get the normalized initial points, lower bound, and upper bound.
            
            % extract the values from the cell array
            x0 = var.x0;
            lb = var.lb;
            ub = var.ub;
            scale = var.scale;
            
            % scale and assign the results to vectors
            [x0_scale, lb_scale, ub_scale] =  SolverUtils.get_var_scale(x0, lb, ub, scale);
        end

        function [name, idx, x, is_bound] = get_param_opt(x_scale, var, lb_scale, ub_scale, tol_bound)
            % Get the parameters for an optimization variable.
            
            idx = var.idx;
            name = var.name;
            lb = var.lb;
            ub = var.ub;
            scale = var.scale;

            tol = tol_bound.*(ub_scale-lb_scale);
            lb_scale_tol = lb_scale+tol;
            ub_scale_tol = ub_scale-tol;
            
            is_bound = (x_scale>lb_scale_tol)&(x_scale<ub_scale_tol);
            x = SolverUtils.get_var_unscale(x_scale, lb, ub, scale);
        end
        
        function [name, idx, x, is_bound] = get_param_fix(var, n_rep)
            % Get the parameters for a fixed variable.
            
            idx = var.idx;
            name = var.name;
            x = var.x;
            
            x = repmat(x, n_rep, 1);
            is_bound = NaN(n_rep, 1);
        end
    end
end