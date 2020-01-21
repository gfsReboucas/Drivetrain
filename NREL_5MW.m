classdef NREL_5MW < Drivetrain 
    %NREL_5MW Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = NREL_5MW
            N_st = 3;
            stage = [Gear_Set Gear_Set Gear_Set];
            for idx = 1:N_st
                stage(idx) =  Gear_Set.NREL_5MW(idx);
            end
            
            P_r = 5.0e3; % [kW], Rated power
            n_r = 12.1; % [1/min.], Input speed
            inp_shaft = Shaft;
            
            m_R = 110.0e3;      J_R = 57231535.0;
            m_G = 1900.0;       J_G = 534.116;
            
            obj@Drivetrain(N_st, stage, P_r, n_r, inp_shaft, m_R, J_R, m_G, J_G);
            
        end
        
    end
end

