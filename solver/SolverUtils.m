classdef SolverUtils < handle
    %% error
    methods(Static, Access = public)
        function err = get_norm(err_vec, wgt_vec, norm)
            % Get the norm on a error vector.
            
            if isfinite(norm)
                err_wgt_vec = abs(err_vec).*(wgt_vec.^(1./norm));
                n_elem = sum(wgt_vec, 2);
                err = (sum(err_wgt_vec.^norm, 2)./n_elem).^(1./norm);
            elseif isinf(norm)
                err = max(abs(err_vec), [], 2);
            else
                error('invalid norm')
            end
            
            idx = any(isfinite(err_vec)==false, 2);
            err(idx) = NaN;
        end
        
        function err = get_percentile(err_vec, percentile)
            % Get the percentile on a error vector.
            
            err = quantile(abs(err_vec), percentile, 2);
            
            idx = any(isfinite(err_vec)==false, 2);
            err(idx) = NaN;
        end
    end
    
    %% clamp
    methods(Static, Access = public)
        function [x_unclamp, lb_unclamp, ub_unclamp] = get_var_unclamp(x, clamp_bnd, lb, ub)
            % Unclamp a variable (sine tranformation).
            
            if clamp_bnd==true
                lb_unclamp = -Inf;
                ub_unclamp = +Inf;
                
                if isinf(lb)&&isinf(ub)
                    x_unclamp = x;
                elseif isinf(lb)
                    x_unclamp = sqrt(x-lb);
                elseif isinf(ub)
                    x_unclamp = sqrt(ub-x);
                else
                    x_unclamp = 2.*(x-lb)./(ub-lb)-1;
                    x_unclamp = 2.*pi+asin(x_unclamp);
                end
            else
                x_unclamp = x;
                lb_unclamp = lb;
                ub_unclamp = ub;
            end
            
            if x_unclamp<=lb_unclamp
                x_unclamp = lb_unclamp+eps;
            end
            
            if x_unclamp>=ub_unclamp
                x_unclamp = ub_unclamp-eps;
            end
        end
        
        function x = get_var_clamp(x_unclamp, clamp_bnd, lb, ub)
            % Clamp a variable (sine tranformation).
            
            if clamp_bnd==true
                if isinf(lb)&&isinf(ub)
                    x = x_unclamp;
                elseif isinf(lb)
                    x = lb+x_unclamp.^2;
                elseif isinf(ub)
                    x = ub-x_unclamp.^2;
                else
                    x = (sin(x_unclamp)+1)./2;
                    x = x.*(ub-lb)+lb;
                end
            else
                x = x_unclamp;
            end
        end
    end
    
    %% scale
    methods(Static, Access = public)
        function [x0_scale, lb_scale, ub_scale] = get_var_scale(x0, lb, ub, scale, norm)
            % Normalize a variable.
            
            switch scale
                case 'lin'
                    fct_scale = @(x) x;
                case 'sqrt'
                    fct_scale = @(x) sqrt(x);
                case 'quad'
                    fct_scale = @(x) x.^2;
                case 'log'
                    fct_scale = @(x) log10(x);
                otherwise
                    error('invalid data')
            end
            
            x0_scale = fct_scale(x0);
            lb_scale = fct_scale(lb);
            ub_scale = fct_scale(ub);
            
            if norm==true
                fct = @(x) (x-lb)./(ub-lb);
                
                x0_scale = fct(x0_scale);
                lb_scale = fct(lb_scale);
                ub_scale = fct(ub_scale);
            end
        end
        
        function x = get_var_unscale(x_scale, lb, ub, scale, norm)
            % Denormalize a variable.
            
            switch scale
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

            if norm==true
                lb_scale = fct_scale(lb);
                ub_scale = fct_scale(ub);
                x_scale = lb_scale+x_scale.*(ub_scale-lb_scale);
            end
            
            x = fct_unscale(x_scale);
        end
    end
    
    %% disp
    methods(Static, Access = public)
        function txt = get_disp_vec(vec)
            % Parse a vector to string.
            
            for i=1:length(vec)
                if islogical(vec(i))
                    txt{i} = mat2str(vec(i));
                elseif isnumeric(vec(i))
                    txt{i} = sprintf('%.3e', vec(i));
                else
                    error('invalid data')
                end
            end
            
            if isscalar(vec)
                txt = txt{:};
            else
                txt = sprintf('[%s]', strjoin(txt, ' ; '));
            end
        end
        
        function txt = get_disp_nan(vec)
            % Parse a vector to NaN.
            
            txt = sprintf('NaN(%d, %d)', size(vec, 1), size(vec, 2));
            
        end
    end
end