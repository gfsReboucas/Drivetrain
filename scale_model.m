classdef scale_model
    properties
        gamma_P (1, :) = 1.0;
        gamma   (1, :) = ones(1, 2);
        DT_ref  (1, :) Drivetrain = Drivetrain;
    end
    
    properties(Dependent)
        DT_sca;
    end
    
    properties(SetAccess = private)
        DT_sca_ (1, :) Drivetrain = Drivetrain;
    end
    
    methods
        function obj = scale_model(gmP, gm, drive_ref)
            if(nargin == 0)
                gmP = 1.0;
                gm = ones(1, 1);
                drive_ref = Drivetrain;
            end
            
            mshaft_ref = drive_ref.main_shaft;
            mshaft_sca = Shaft(mshaft_ref.d*obj.gamma(3, :), ...
                               mshaft_ref.L*obj.gamma(4, :));
            
            ref_stage = zeros(1, 3);
            sca_stage = zeros(1, 3);
                        
            for idx = 1:3
                ref_stage(idx) = drive_ref.stage(idx);
                stage_shaft_ref = ref_stage.shaft;
                stage_shaft_sca = Shaft(stage_shaft_ref.d*obj.gamma(3, :), ...
                                        stage_shaft_ref.L*obj.gamma(4, :));
                
                sca_stage(idx) = Gear_Set(ref_stage(idx).configuration, ...
                                          ref_stage(idx).m_n*obj.gamma(1, :), ...
                                          ref_stage(idx).alpha_n, ...
                                          ref_stage(idx).z, ...
                                          ref_stage(idx).b*obj.gamma(2, :), ...
                                          ref_stage(idx).x, ...
                                          ref_stage(idx).beta, ...
                                          ref_stage(idx).k, ...
                                          ref_stage(idx).bore_ratio, ...
                                          ref_stage(idx).N_p, ...
                                          ref_stage(idx).a_w*obj.gamma(1, :), ...
                                          ref_stage(idx).type, ...
                                          ref_stage(idx).bearing, ...
                                          ref_stage(idx).N_p, ...
                                          stage_shaft_sca);
            end
            
        end
        
    end
    
    %% Get methods:
    methods
        function val = get.DT_sca(gm)
            
        end
        
    end
    
    
end
