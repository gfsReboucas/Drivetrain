classdef scaling_factor
    properties
        name;
        value;
    end
    
    methods
        function obj = scaling_factor(par, val)
            if(nargin == 0)
                par = "*";
                val = NaN;
            end
            
            if(numel(par) ~= numel(val))
                error("Name and value should have the same number of elements.");
            end
            
            if(isrow(val))
                val = val';
            end
            
            obj.name = par;
            obj.value = val;
        end
        
        function disp(obj)
            if(isnumeric([obj.name{:}]))
                formatSpec = "%d:\t%s\n";
            else
                formatSpec = "%s:\t%.3e\n";
            end
            
            for idx = 1:numel(obj.name)
                fprintf(formatSpec, obj.name{idx}, obj.value(idx));
            end
        end
        
        function val = length(obj)
            val = length(obj.value);
        end
        
%         function val = numel(obj) 
%         % removed because it was messing up subsref
%             val = numel(obj.value);
%         end
        
        function val = size(obj)
            val = [length(obj) 1];
        end
        
    end
    
    %% Set methods:
    methods
        function obj = subsasgn(obj, par, varargin)
            if isequal(obj,[])
                % obj = ClassName.empty;
            end
            
            if(isrow(varargin{:}))
                varargin{:} = varargin{:}';
            end
            
            switch(par(1).type)
                case '.'
                    if length(par) == 1
                        obj.(par.subs) = varargin{:};
                    elseif((length(par) == 2) && (strcmp(par(2).type, '()')))
                        obj.(par(1).subs)(par(2).subs{:}) = varargin{:};
                    else
                        obj = builtin('subsasgn', obj, par, varargin{:});
                    end
                case '()'
                    par_idx = par.subs{:};
                    par_idx(par_idx == "*") = [];

                    if(numel(par_idx) ~= numel(varargin{:}))
                        error("Name and value should have the same number of elements.");
                    end

                    if(isnumeric(par_idx))
                        if((min(par_idx) >= 1) && (max(par_idx) <= numel(obj)))
                            idx = par_idx;
                        else
                            error("Indexes out of range.");
                        end
                    elseif(isstring(par_idx))
                        idx = ismember(obj.name, par_idx);
                    end
                    
                    obj.value(idx) = varargin{:}';
                otherwise
                    error("Operator [%s] not implemented yet.", par(1).type);
            end
        end
    end
    
    %% Get methods:
    methods
        function val = get.value(obj)
            % based on: (10/02/2020)
            % https://www.mathworks.com/matlabcentral/answers/337824-why-jsonencode-with-containers-map-doesn-t-preserve-order-of-key-value-pairs
            val = obj.value;
        end
        
        function val = get.name(obj)
            val = obj.name;
        end
        
        function varargout = subsref(obj, par)
            
            switch(par(1).type)
                case '.'
                    if(length(par) == 1)
                        val = obj.(par.subs);
                    elseif((length(par) == 2) && (strcmp(par(2).type,'()')))
                        val = obj.(par(1).subs)(par(2).subs{:});
                    else
                        val = builtin('subsref', obj, par);
                    end
                case '()'
                    par_idx = par.subs{:};

                    if(isnumeric(par_idx))
                        if((min(par_idx) >= 1) && (max(par_idx) <= length(obj)))
                            idx = par_idx;
                        else
                            error("Indexes out of range.");
                        end
                    elseif(islogical(par_idx))
                        if((numel(par_idx) > 1) && (numel(par_idx) <= length(obj)))
                            idx = par_idx;
                        else
                            error("Indexes out of range.");
                        end
                    elseif(isstring(par_idx))
                        idx = ismember(obj.name, par_idx);
                    end

                    val = obj.value(idx);
                otherwise
                    error("Operator [%s] not implemented yet.", par(1).type);
            end
            
            [varargout{1:nargout}] = val;
        end
        
        function ind = end(obj, k, n)
            szd = size(obj);
            if k < n
                ind = szd(k);
            else
                ind = prod(szd(k:end));
            end
        end
    end
    
    %% Static method:
    methods(Static)
        function gamma = generate_gamma(n, N)
            %GENERATE_GAMMA generates an array of values from power(10, -n) to 1 with N
            % elements.
            %
            
            exp_array = -n:0;
            specific = 10.^exp_array;
            
            p = ceil(N/n);
            gm_lin = zeros(n, p);
            gm_log = zeros(n, p);
            
            for idx = 1:n
                gm_lin(idx, :) = linspace(specific(idx) , specific(idx + 1) , p);
                gm_log(idx, :) = logspace(exp_array(idx), exp_array(idx + 1), p);
            end
            
            gm_lin = reshape(gm_lin, 1, n*p);
            gm_log = reshape(gm_log, 1, n*p);
            
            gm = sort([gm_lin, gm_log]);
            gamma = unique(gm(1:2:end));
        end
    end
    
end
