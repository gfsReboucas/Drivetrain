function out = process_varargin(vargin, standard)
    vargin = struct(vargin{:});
    standard = struct(standard{:});
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
    
end
