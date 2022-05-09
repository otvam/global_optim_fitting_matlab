classdef SolverUtils < handle
    % Static class with various utils.
    %
    %    Computing error metrics.
    %    Scaling of variables.
    %
    %    Thomas Guillod.
    %    2021-2022 - BSD License.
    
    %% error
    methods(Static, Access = public)
        function out = get_norm(err_mat, wgt_mat, norm)
            % Get the norm on an error matrix with specified weights.
            
            % check
            idx_nan = SolverUtils.get_err_data(err_mat, wgt_mat);
            
            % compute the norm
            if isfinite(norm)
                err_wgt_mat = abs(err_mat).*(wgt_mat.^(1./norm));
                n_elem = sum(wgt_mat, 1);
                out = (sum(err_wgt_mat.^norm, 1)./n_elem).^(1./norm);
            elseif isinf(norm)
                out = max(abs(err_mat), [], 1);
            else
                error('invalid norm')
            end
            
            % replace invalid data
            out(idx_nan) = NaN;
        end
        
        function out = get_error(err_mat, wgt_mat, type)
            % Get simple weighted error metrics for an error matrix.
            
            % check
            idx_nan = SolverUtils.get_err_data(err_mat, wgt_mat);
            
            % compute the value
            switch type
                case 'avg'
                    out = sum(err_mat.*wgt_mat, 1)./sum(wgt_mat, 1);
                case 'min'
                    out = min(err_mat, [], 1);
                case 'max'
                    out = max(err_mat, [], 1);
                otherwise
                    error('invalid data')
            end
            
            % replace invalid data
            out(idx_nan) = NaN;
        end
        
        function out = get_percentile(err_mat, wgt_mat, percentile)
            % Get simple the weighted percentile for an error matrix.
            
            % check
            idx_nan = SolverUtils.get_err_data(err_mat, wgt_mat);
            
            % compute the weighted percentile
            for i=1:size(err_mat, 2)
                % extract
                err_vec = err_mat(:,i);
                wgt_vec = wgt_mat(:,i);
                
                % remove duplicated values
                [err_vec, idx, idx_rev] = unique(err_vec);
                wgt_vec = accumarray(idx_rev, wgt_vec, [], @sum);
                tmp_vec = 1:length(idx);
                
                % get the cumulative weights between zero and one
                out_vec = (cumsum(wgt_vec)-0.5.*wgt_vec)./sum(wgt_vec);
                
                % compute the weighted percentile
                if percentile<=min(out_vec)
                    err = min(err_vec);
                elseif percentile>=max(out_vec)
                    err = max(err_vec);
                else
                    idx = interp1(out_vec, tmp_vec, percentile);
                    err = interp1(tmp_vec, err_vec, idx);
                end
                
                % assign
                out(i) = err;
            end
            
            % replace invalid data
            out(idx_nan) = NaN;
        end
        
        function idx_nan = get_err_data(err_mat, wgt_mat)
            % Check the validity of an error matrix with specified weights.
            
            assert(all(size(err_mat)==size(wgt_mat)), 'invalid data')
            idx_nan = any(isfinite(err_mat)==false, 1)|any(isfinite(wgt_mat)==false, 1);
        end
    end
    
    %% scaling
    methods(Static, Access = public)
        function [x_unclamp, lb_unclamp, ub_unclamp] = get_var_unclamp(x, lb, ub, clamp_bnd)
            % Transform bounded a variable into a unconstrained variable with sine transformation.
            
            if clamp_bnd==true
                % the new variable is unconstrained
                lb_unclamp = -Inf;
                ub_unclamp = +Inf;
                
                if isinf(lb)&&isinf(ub)
                    % do nothing as the variable is already unconstrained
                    x_unclamp = x;
                elseif isinf(lb)
                    % quadratic transformation for single-sided bound
                    x_unclamp = sqrt(x-lb);
                elseif isinf(ub)
                    % quadratic transformation for single-sided bound
                    x_unclamp = sqrt(ub-x);
                else
                    % sine transformation for double-sided bounds
                    x_unclamp = 2.*(x-lb)./(ub-lb)-1;
                    x_unclamp = asin(x_unclamp);
                    
                    % shift to avoid numerical issue around zero
                    x_unclamp = 2.*pi+x_unclamp;
                end
            else
                x_unclamp = x;
                lb_unclamp = lb;
                ub_unclamp = ub;
            end
            
            % avoid numerical issue at the bounds
            if x_unclamp<=lb_unclamp
                x_unclamp = lb_unclamp+eps;
            end
            if x_unclamp>=ub_unclamp
                x_unclamp = ub_unclamp-eps;
            end
        end
        
        function x = get_var_clamp(x_unclamp, lb, ub, clamp_bnd)
            % Transform a unconstrained variable into a bounded variable with sine transformation.
            
            if clamp_bnd==true
                if isinf(lb)&&isinf(ub)
                    % do nothing as the variable is already unconstrained
                    x = x_unclamp;
                elseif isinf(lb)
                    % quadratic transformation for single-sided bound
                    x = lb+x_unclamp.^2;
                elseif isinf(ub)
                    % quadratic transformation for single-sided bound
                    x = ub-x_unclamp.^2;
                else
                    % sine transformation for double-sided bounds
                    x = (sin(x_unclamp)+1)./2;
                    x = x.*(ub-lb)+lb;
                end
            else
                x = x_unclamp;
            end
        end
        
        function x_trf = get_var_trf(x, trf, is_revert)
            % Variable transformation (or reverse transformation).
            
            % get the scaling functions
            switch trf
                case 'lin'
                    fct_scale = @(x) x;
                    fct_unscale = @(x) x;
                case 'sqrt'
                    fct_scale = @(x) sqrt(x);
                    fct_unscale = @(x) x.^2;
                case 'quad'
                    fct_scale = @(x) x.^2;
                    fct_unscale = @(x) sqrt(x);
                case 'log'
                    fct_scale = @(x) log10(x);
                    fct_unscale = @(x) 10.^x;
                otherwise
                    error('invalid data')
            end
            
            % transform the variable
            if is_revert==true
                x_trf = fct_unscale(x);
            else
                x_trf = fct_scale(x);
            end
        end
        
        function x_norm = get_var_norm(x, lb, ub, norm, is_revert)
            % Variable normalization (or denormalization).
            
            if norm==true
                if is_revert==true
                    x_norm = lb+x.*(ub-lb);
                else
                    x_norm = (x-lb)./(ub-lb);
                end
            end
        end
    end
end