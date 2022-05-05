classdef SolverCache < handle
    %% properties
    properties (SetAccess = private, GetAccess = private)
        fct_opt
        use_cache
        vec_cache
        n_cache
        tol_cache
        
        i_cache
        x_val
        err_val
        wgt_val
    end
    
    %% public
    methods (Access = public)
        function self = SolverCache(fct_opt, use_cache, vec_cache, n_cache, tol_cache)
            % set data
            self.fct_opt = fct_opt;
            self.use_cache = use_cache;
            self.vec_cache = vec_cache;
            self.n_cache = n_cache;
            self.tol_cache = tol_cache;
            
            % init data
            self.get_clear();
        end
        
        function [err, wgt] = get_eval(self, x)
            if self.use_cache==true
                [err, wgt] = self.get_eval_cache(x);
            else
                [err, wgt] = self.get_eval_fct(x);
            end
        end
        
        function get_clear(self)
            self.i_cache = 0;
            self.x_val = [];
            self.err_val = [];
            self.wgt_val = [];
        end
    end
    
    %% private api
    methods (Access = private)
        function [err, wgt] = get_eval_fct(self, x)
            % Evaluate the function.
            
            if self.vec_cache
                [err, wgt] = self.fct_opt(x);
            else
                for i=1:size(x, 1)
                    [err(i,:), wgt(i,:)] = self.fct_opt(x(i,:));
                end
            end
        end
        
        function [err, wgt] = get_eval_cache(self, x)
            % Cache call.
            
            for i=1:size(x, 1)
                idx_select(i) = self.get_find(x(i,:));
            end
            [err, wgt] = self.get_value(x, idx_select);
            self.get_update(x, err, wgt, idx_select);
        end
        
        function idx_select = get_find(self, x)
            % Find if an element is in the cache.
            
            if self.i_cache==0
                idx_select = NaN;
            else
                idx = isnan(x)&isnan(self.x_val);
                err = x-self.x_val;
                err(idx) = 0;
                err = vecnorm(err, 2, 2);
                
                idx_remove = find(err<self.tol_cache);
                
                if isempty(idx_remove)
                    idx_select = NaN;
                else
                    [~, idx_select] = min(err(idx_remove));
                    idx_select = idx_remove(idx_select);
                end
            end
        end
        
        function [err, wgt] = get_value(self, x, idx_select)
            % Get the values.
            
            % get from evaluation
            idx_compute = isnan(idx_select);
            if nnz(idx_compute)==0
                err_compute = [];
                wgt_compute = [];
            else
                [err_compute, wgt_compute] = self.get_eval_fct(x(idx_compute,:));
            end
            
            % get from cache
            idx_cache = isfinite(idx_select);
            idx_access = idx_select(idx_cache);
            err_cache = self.err_val(idx_access,:);
            wgt_cache = self.wgt_val(idx_access,:);

            % assemble cache
            err(idx_compute,:) = err_compute;
            wgt(idx_compute,:) = wgt_compute;
            err(idx_cache,:) = err_cache;
            wgt(idx_cache,:) = wgt_cache;
        end
        
        function get_update(self, x, err, wgt, idx_select)
            % Update the cache.
                        
            % remove points
            idx_remove = isfinite(idx_select);
            idx_remove = idx_select(idx_remove);
            self.x_val(idx_remove,:) = [];
            self.err_val(idx_remove,:) = [];
            self.wgt_val(idx_remove,:) = [];
            self.i_cache = self.i_cache-length(idx_remove);
            
            % add points
            self.x_val = [self.x_val ; x];
            self.err_val = [self.err_val ; err];
            self.wgt_val = [self.wgt_val ; wgt];
            self.i_cache = self.i_cache+length(idx_select);
            
            % limit cache size
            if self.i_cache>self.n_cache
                idx = 1:(self.i_cache-self.n_cache);
                self.i_cache = self.n_cache;
                
                self.x_val(idx,:) = [];
                self.err_val(idx,:) = [];
                self.wgt_val(idx,:) = [];
            end
        end
    end
end