function out = process_varargin(vargin, standard)
    standard = struct(standard{:});
    
    if(~isempty(vargin))
        vargin = struct(vargin{:});
        out = standard;
        
        field_v = fields(vargin);
        field_s = fields(standard);
        
        lv = length(field_v);
        ls = length(field_s);
        
        if(lv > ls)
            warning('varargin has more parameters than standard');
        end
        
        for idx = 1:lv
            f = field_v{idx};
            if(isfield(standard, f))
                out.(f) = vargin.(f);
            end
        end
    else
        out = standard;
    end
    
end
