classdef Twin_Designer < Drivetrain
    properties
        gamma_P;
        gamma_n;
        gamma scaling_factor;
    end
    
    methods
        function obj = Twin_Designer(name, gm_P, gm_n)
            if(strcmp(name, 'NREL_5MW'))
                obj = NREL_5MW();
            elseif(strcmp(name, 'DTU_10MW'))
                obj = DTU_10MW();
            else
                error('name = [%s] not recognized.', upper(name));
            end
            
            obj.gamma_P = gm_P;
            obj.gamma_n = gm_n;
        end
    end
end
