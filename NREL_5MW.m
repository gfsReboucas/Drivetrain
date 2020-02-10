classdef NREL_5MW < Drivetrain % matlab.mixin.SetGet &
    %NREL_5MW This class contains some of the properties of the NREL 5MW
    % wind turbine gearbox proposed by Nejad et. al. [1]. More information
    % regarding the NREL 5MW wind turbine can be found at [2].
    %
    % [1] Nejad, A. R., Guo, Y., Gao, Z., Moan, T. (2016). Development of a
    % 5 MW reference gearbox for offshore wind turbines. Wind Energy. 
    % https://doi.org/10.1002/we.1884
    % [2] Jonkman, J., Butterfield, S., Musial, W., & Scott, G. Definition
    % of a 5-MW Reference Wind Turbine for Offshore System Development. 
    % doi:10.2172/947422.
    % [3] Anaya-Lara, O., Tande, J.O., Uhlen, K., Merz, K. and Nejad, A.R.
    % (2019). Modelling and Analysis of Drivetrains in Offshore Wind
    % Turbines. In Offshore Wind Energy Technology (eds O. Anaya-Lara, J.O.
    % Tande, K. Uhlen and K. Merz). doi:10.1002/9781119097808.ch3
    %
    % written by:
    % Geraldo Rebouças
    % - Geraldo.Reboucas@ntnu.no OR
    % - gfs.reboucas@gmail.com
    %
    % Postdoctoral Fellow at:
    %    Norwegian University of Science and Technology, NTNU
    %    Department of Marine Technology, IMT
    %    Marine System Dynamics and Vibration Lab, MD Lab
    %    https://www.ntnu.edu/imt/lab/md-lab
    %    
    % see also DRIVETRAIN.
    
    properties
        gamma containers.Map;
    end
    
    properties(Dependent)
        % Stage 1:
        m_n1; % Normal module, [mm]
        b_1;  % Face width, [mm]
        d_1;  % Diameter of the output shaft, [mm]
        L_1;  % Length of the output shaft, [mm]
        % Stage 2:
        m_n2; % Normal module, [mm]
        b_2;  % Face width, [mm]
        d_2;  % Diameter of the output shaft, [mm]
        L_2;  % Length of the output shaft, [mm]
        % Stage 3:
        m_n3; % Normal module, [mm]
        b_3;  % Face width, [mm]
        d_3;  % Diameter of the output shaft, [mm]
        L_3;  % Length of the output shaft, [mm]
        
        J_R;  % Mass moment of inertia of the rotor, [kg-m^2] 
        J_G;  % Mass moment of inertia of the rotor, [kg-m^2]
        
        d_s;  % Diameter of the main shaft, [mm]
        L_s;  % Length of the main shaft, [mm]
    end

    methods
        function obj = NREL_5MW
            N_st = 3;
            stage = [Gear_Set Gear_Set Gear_Set];
            
            for idx = 1:N_st
                stage(idx) =  NREL_5MW.gear_stage(idx);
            end
            
            P_r = 5.0e3; % [kW], Rated power
            n_r = 12.1; % [1/min.], Input speed
            inp_shaft = NREL_5MW.shaft(0);
            
            m_R = 110.0e3;      J_R = 57231535.0; % according to [1, 3]
            m_G = 1900.0;       J_G = 534.116;
            
%             r_R: 3.226e+04	h_R: 4.297e+00 [mm]
%             r_G: 7.498e+02	h_G: 1.374e+02 [mm]
            
            obj@Drivetrain(N_st, stage, P_r, n_r, inp_shaft, m_R, J_R, m_G, J_G);
            obj.dynamic_model =  "Kahraman_1994";
            
            %          Normal module, Face width, Length, Diameter (output shaft)
            key_set = ["m_n1"       , "b_1"     , "d_1" , "L_1", ... % Stage 01
                       "m_n2"       , "b_2"     , "d_2" , "L_2", ... % Stage 02
                       "m_n3"       , "b_3"     , "d_3" , "L_3", ... % Stage 03
                       "J_R",                                    ... %   M. M. Inertia (rotor)
                       "J_G",                                    ... %   M. M. Inertia (generator)
                                                  "d_s" , "L_s"];    % Main shaft
                   
            val_set = ones(size(key_set));
            
            obj.gamma = containers.Map(key_set, val_set);
        end
    end
    
    methods(Static)
        function brg = bearing(idx)
            %BEARING returns the bearing sets of each stage of the NREL 
            % 5 MW wind turbine drivetrain according to [1], tables VI and
            % IX.
            %
            
            switch(idx)
                case 0 % Main shaft
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    INP_A     = Bearing("INP_A"    , "CARB", 0.0   , 1.50e10, 1.50e10, 0.0, 5.0e6, 5.0e6, 1750.0, 1250.0, 375.0);
                    INP_B     = Bearing("INP_B"    , "SRB" , 4.06e8, 1.54e10, 1.54e10, 0.0, 0.0  , 0.0  , 1220.0,  750.0, 365.0);
                    
                    brg = [INP_A, INP_B];
                    
                case 1
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    PL_A      = Bearing("PL_A"     , "CRB",  9.1e4,  9.4e9,   3.2e9,   0.0, 1.4e6, 4.5e6, 600.0,  400.0,  272.0);
                    PL_B      = Bearing("PL_B"     , "CRB",  9.1e4,  9.4e9,   3.2e9,   0.0, 1.4e6, 4.5e6, 600.0,  400.0,  272.0);
                    PLC_A     = Bearing("PLC_A"    , "SRB",  6.6e4,  1.7e9,   1.1e9,   0.0, 5.6e5, 1.3e5, 1030.0,  710.0, 315.0);
                    PLC_B     = Bearing("PLC_B"    , "CRB",  6.6e7,  1.7e9,   1.1e9,   0.0, 5.6e5, 1.3e5, 1220.0, 1000.0, 128.0);
                    
                    brg = [PL_A,  PL_B, ... % Planet
                           PLC_A, PLC_B];   % Carrier
                    
                case 2
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    IMS_PL_A  = Bearing("IMS_PL_A" , "CRB",  9.1e4,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 520.0,  380.0,  140.0);
                    IMS_PL_B  = Bearing("IMS_PL_B" , "CRB",  9.1e4,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 520.0,  380.0,  140.0);
                    IMS_PLC_A = Bearing("IMS_PLC_A", "CARB", 9.1e4,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 1030.0, 710.0,  236.0);
                    IMS_PLC_B = Bearing("IMS_PLC_B", "CRB" , 9.1e7,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4,  870.0, 600.0,  155.0);
                    
                    brg = [IMS_PL_A,  IMS_PL_B, ... % Planet
                           IMS_PLC_A, IMS_PLC_B];   % Carrier
                    
                case 3
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    IMS_A     = Bearing("IMS_A"    , "CRB",  0.0,    6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 360.0,  200.0,   98.0);
                    IMS_B     = Bearing("IMS_B"    , "TRB",  7.4e7,  5.0e8,   5.0e8,   0.0, 1.6e6, 1.8e6, 460.0,  200.0,  100.0);
                    IMS_C     = Bearing("IMS_C"    , "TRB",  7.8e7,  7.4e8,   3.3e8,   0.0, 1.1e6, 2.5e6, 460.0,  200.0,  100.0);
                    HS_A      = Bearing("HS_A"     , "CRB",  1.3e8,  8.2e8,   8.2e8,   0.0, 1.7e5, 1.0e6, 500.0,  400.0,  100.0);
                    HS_B      = Bearing("HS_B"     , "TRB",  6.7e7,  8.0e8,   1.3e8,   0.0, 1.7e5, 1.0e6, 550.0,  410.0,   86.0);
                    HS_C      = Bearing("HS_C"     , "TRB",  8.0e7,  1.0e9,   7.3e7,   0.0, 1.7e5, 1.0e6, 550.0,  410.0,   86.0);
                    
                    brg = [HS_A,  HS_B,  HS_C, ... % Pinion
                           IMS_A, IMS_B, IMS_C];   % Wheel

                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", idx);
            end
        end
        
        function sft = shaft(idx)
            %SHAFT returns the shafts of each stage of the NREL 5 MW
            % wind turbine drivetrain according to [1]. Based mainly on the
            % SIMPACK simulation provided by its first author.
            %
            
            switch(idx)
                case 0 % LSS
%                     d_0 = 1000.0;
%                     L_0 = 3000.0;
                    d_0 = 700.0;
                    L_0 = 2000.0;

                    sft = Shaft(d_0, L_0);
                    
                case 1 % ISS
%                     d_1 = 800.0;
%                     L_1 = 750.0;
                    d_1 = 533.0;
                    L_1 = 500.0;

                    sft = Shaft(d_1, L_1);
                    
                case 2 % HS-IS
%                     d_2 = 500.0;
%                     L_2 = 1000.0;
                    d_2 = 333.0;
                    L_2 = 666.0;

                    sft = Shaft(d_2, L_2);
                    
                case 3 % HSS
%                     d_3 = 500.0;
%                     L_3 = 1500.0;
                    d_3 = 333.0;
                    L_3 = 1000.0;

                    sft = Shaft(d_3, L_3);
                    
                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", idx);
            end
        end
        
        function g_set = gear_stage(idx)
            %GEAR_STAGE returns the gear stages of the NREL 5 MW wind
            % turbine drivetrain according to [4], Table V. The values for
            % the tip alteration coefficients were taken from KISSsoft.
            %

            alpha_n = 20.0;        % [deg.],   Pressure angle (at reference cylinder)
            rack_type = "A";       % [-],      Type of the basic rack from A to D

            switch(idx)
                case 1
                    p_1    = 3;            % [-],    Number of planet gears
                    m_n1   = 45.0;         % [mm],   Normal module
                    beta_1 = 0.0;          % [deg.], Helix angle (at reference cylinder)
                    b_1    = 491.0;        % [mm],   Face width
                    a_w1   = 863.0;        % [mm],   Center distance
                    z_s1   =  19;          % [-],    Number of teeth (sun)    [WHEEL]
                    z_p1   =  17;          % [-],    Number of teeth (planet) [PINION]
                    z_r1   =  56;          % [-],    Number of teeth (ring)
                    x_s1   =  0.617;       % [-],    Profile shift coefficient (sun)
                    x_p1   =  0.802;       % [-],    Profile shift coefficient (planet)
                    x_r1   = -0.501;       % [-],    Profile shift coefficient (ring)
                    k_s1   = -10.861/m_n1; % [-],    Tip alteration coefficient (sun)
                    k_p1   = -10.861/m_n1; % [-],    Tip alteration coefficient (planet)
                    k_r1   =   0.0;        % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs1 = 80.0/171.0;
                    bore_Rp1 = 80.0/153.0;
                    bore_Rr1 = 1.2;
                    
                    z_1 = [z_s1 z_p1 z_r1];
                    x_1 = [x_s1 x_p1 x_r1];
                    k_1 = [k_s1 k_p1 k_r1];
                    bore_R1 = [bore_Rs1 bore_Rp1 bore_Rr1];
%                                 Sun,   Planet,   Ring,   Carrier
                    bearing_1 = NREL_5MW.bearing(1);
                    shaft_1 = NREL_5MW.shaft(1);
                    
                    g_set = Gear_Set("planetary", m_n1, alpha_n, z_1, b_1, x_1, beta_1, k_1, bore_R1, p_1, a_w1, rack_type, bearing_1, shaft_1);

                case 2
                    p_2    = 3;          % [-],    Number of planet gears
                    m_n2   = 21.0;       % [mm],   Normal module
                    beta_2 = 0.0;        % [deg.], Helix angle (at reference cylinder)
                    b_2    = 550.0;      % [mm],   Face width
                    a_w2   = 584.0;      % [mm],   Center distance
                    z_s2   =  18;        % [-],    Number of teeth (sun)    [PINION]
                    z_p2   =  36;        % [-],    Number of teeth (planet) [WHEEL]
                    z_r2   =  93;        % [-],    Number of teeth (ring)
                    x_s2   = 0.389;      % [-],    Profile shift coefficient (sun)
                    x_p2   = 0.504;      % [-],    Profile shift coefficient (planet)
                    x_r2   = 0.117;      % [-],    Profile shift coefficient (ring)
                    k_s2   = -1.75/m_n2; % [-],    Tip alteration coefficient (sun)
                    k_p2   = -1.75/m_n2; % [-],    Tip alteration coefficient (planet)
                    k_r2   =   0.0;      % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs2 = 100.0/189.0;
                    bore_Rp2 = 95.0/189.0;
                    bore_Rr2 = 1.2;

                    z_2 = [z_s2 z_p2 z_r2];
                    x_2 = [x_s2 x_p2 x_r2];
                    k_2 = [k_s2 k_p2 k_r2];
                    bore_R2 = [bore_Rs2 bore_Rp2 bore_Rr2];
%                                 Sun,   Planet,   Ring,   Carrier
                    bearing_2 = NREL_5MW.bearing(2);
                    shaft_2 = NREL_5MW.shaft(2);
                    
                    g_set = Gear_Set("planetary", m_n2, alpha_n, z_2, b_2, x_2, beta_2, k_2, bore_R2, p_2, a_w2, rack_type, bearing_2, shaft_2);

                case 3
                    p_3    = 1;            % [-],    Number of planet gears
                    m_n3   = 14.0;         % [mm],   Normal module
                    beta_3 = 10.0;         % [deg.], Helix angle (at reference cylinder)
                    b_3    = 360.0;        % [mm],   Face width
                    a_w3   = 861.0;        % [mm],   Center distance
                    z_13   =  24;          % [-],    Number of teeth (pinion)
                    z_23   =  95;          % [-],    Number of teeth (wheel)
                    x_13   =  0.480;       % [-],    Profile shift coefficient (pinion)
                    x_23   =  0.669;       % [-],    Profile shift coefficient (wheel)
                    k_13   =  -0.938/14.0; % [-],    Tip alteration coefficient (pinion)
                    k_23   =  -0.938/14.0; % [-],    Tip alteration coefficient (wheel)
                    
                    bore_R13 = 1809.0/3086.0;
                    bore_R23 = 3385.0/9143.0;

                    z_3 = [z_13 z_23];
                    x_3 = [x_13 x_23];
                    k_3 = [k_13 k_23];
                    bore_R3 = [bore_R13 bore_R23];
%                                 Pinion, Wheel
                    bearing_3 = NREL_5MW.bearing(3);
                    shaft_3 = NREL_5MW.shaft(3);
                    
                    g_set = Gear_Set("parallel", m_n3, alpha_n, z_3, b_3, x_3, beta_3, k_3, bore_R3, p_3, a_w3, rack_type, bearing_3, shaft_3);
                    
                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", idx);
            end
        end
    end
    
    %% Set methods:
    methods
%         function obj = set.gamma(obj, val)
%             obj.gamma = containers.Map(obj.gamma.keys, val);
%         end
        
        function set.m_n1(obj, val)
            obj.stage(1).m_n = val;
        end
        
        function set.b_1(obj, val)
            obj.stage(1).b = val;
        end
        
        function set.d_1(obj, val)
            obj.stage(1).out_shaft.d = val;
        end
        
        function set.L_1(obj, val)
            obj.stage(1).out_shaft.L = val;
        end
        
        function set.m_n2(obj, val)
            obj.stage(2).m_n = val;
        end
        
        function set.b_2(obj, val)
            obj.stage(2).b = val;
        end
        
        function set.d_2(obj, val)
            obj.stage(2).out_shaft.d = val;
        end
        
        function set.L_2(obj, val)
            obj.stage(2).out_shaft.L = val;
        end
        
        function set.m_n3(obj, val)
            obj.stage(3).m_n = val;
        end
        
        function set.b_3(obj, val)
            obj.stage(3).b = val;
        end
        
        function set.d_3(obj, val)
            obj.stage(3).out_shaft.d = val;
        end
        
        function set.L_3(obj, val)
            obj.stage(3).out_shaft.L = val;
        end
        
        function set.J_R(obj, val)
            obj.J_Rotor = val;
        end
        
        function set.J_G(obj, val)
            obj.J_Gen = val;
        end
        
        function set.d_s(obj, val)
            obj.main_shaft.d = val;
        end
        
        function set.L_s(obj, val)
            obj.main_shaft.L = val;
        end
        
    end
    %% Get methods:
    methods
%         function val = get.gamma(obj)
%             val = obj.gamma;
%         end
        
        function val = get.m_n1(obj)
            val = obj.stage(1).m;
        end
        
        function val = get.b_1(obj)
            val = obj.stage(1).b;
        end
        
        function val = get.d_1(obj)
            val = obj.stage(1).out_shaft.d;
        end
        
        function val = get.L_1(obj)
            val = obj.stage(1).out_shaft.L;
        end
        
        function val = get.m_n2(obj)
            val = obj.stage(2).m;
        end
        
        function val = get.b_2(obj)
            val = obj.stage(2).b;
        end
        
        function val = get.d_2(obj)
            val = obj.stage(2).out_shaft.d;
        end
        
        function val = get.L_2(obj)
            val = obj.stage(2).out_shaft.L;
        end
        
        function val = get.m_n3(obj)
            val = obj.stage(3).m_n;
        end
        
        function val = get.b_3(obj)
            val = obj.stage(3).b;
        end
        
        function val = get.d_3(obj)
            val = obj.stage(3).out_shaft.d;
        end
        
        function val = get.L_3(obj)
            val = obj.stage(3).out_shaft.L;
        end
        
        function val = get.J_R(obj)
            val = obj.J_Rotor;
        end
        
        function val = get.J_G(obj)
            val = obj.J_Gen;
        end
        
        function val = get.d_s(obj)
            val = obj.main_shaft.d;
        end
        
        function val = get.L_s(obj)
            val = obj.main_shaft.L;
        end
        
    end
    
    %% Scaling calculations:
    methods
        function obj_sca = scale_all(obj_ref, gamma_P, gamma_n, gamma)
            %SCALE_ALL returns a NREL_5MW object scaled by the
            % factors gamma_P and gamma_n for its rated power and rotor
            % speed and gamma for normal module, face width, shaft
            % dimensions (diameter and length), mass and mass moment of
            % inertia of rotor and generator.
            %
            % current limitations/constraints:
            % - shaft diameter is proportional to torque;
            % - dimensions of the main shaft are proportional to the
            % output shaft in stage 03;
            %
            
            if((length(gamma) ~= 18) || (~isa(gamma, "containers.Map")))
                error("gamma must be a container with 18 elements.");
            end
            
            key_set = obj_ref.gamma.keys;
            
            gamma_mR  = gamma("m_R");
            gamma_JR  = gamma("J_R");
            gamma_mG  = gamma("m_G");
            gamma_JG  = gamma("J_G");
%             gamma_Ls  = gamma("L_s");
%             gamma_ds  = gamma("d_s");
            
            gamma_Torque = gamma_P/gamma_n; % Applied torque
            
            % Main shaft:
%             gamma_Ls = gamma("L_3"); % length
            gamma_ds = nthroot(gamma_Torque, 3.0);
            
            % Shaft diameters are scaled in the same way
            idx_dd = find(contains(key_set, "d"));
            
            for idx = idx_dd
                gamma(key_set(idx)) = gamma_ds;
            end
            
            stage_sca = [Gear_Set, Gear_Set, Gear_Set];
            
            idx_LL = find(contains(key_set, "L"));
            
            gamma_dd = cell2mat(values(gamma, {key_set{idx_dd}}))';%#ok<*FNDSB>
            gamma_LL = cell2mat(values(gamma, {key_set{idx_LL}}))';
            
            for idx = 1:obj_ref.N_stg
                jdx = 4*idx + (-3:0);
                gamma_stage = cell2mat(values(gamma, {key_set{jdx}}));
                
                stage_sca(idx) = obj_ref.stage(idx).scale_aspect(gamma_stage, "Gear_Set");
            end
            
            idx = idx + 1;
            
            main_shaft_sca = Shaft(obj_ref.main_shaft.d*gamma_dd(idx), ...
                                   obj_ref.main_shaft.L*gamma_LL(idx));
            
            obj_sca = Drivetrain(obj_ref.N_stg, ...
                                 stage_sca, ...
                                 obj_ref.P_rated *gamma_P, ...
                                 obj_ref.n_rotor *gamma_n, ...
                                 main_shaft_sca, ...
                                 obj_ref.m_Rotor*gamma_mR, ...
                                 obj_ref.J_Rotor*gamma_JR, ...
                                 obj_ref.m_Gen  *gamma_mG, ...
                                 obj_ref.J_Gen  *gamma_JG);
                             
            obj_sca.dynamic_model = obj_ref.dynamic_model;
        end
        
        function obj_sca = scale_aspect(obj_ref, gamma_P, gamma_n, gamma, aspect)
            %SCALE_ASPECT returns a NREL_5MW object scaled by the factors
            % gamma_P and gamma_n for its rated power and rotor speed. The
            % argument gamma can be used to scale the NREL_5MW
            % considering various aspects.
            %
            
            %          Normal module, Face width, Length, Diameter (output shaft)
            key_set = ["m_n1"       , "b_1"     , "d_1" , "L_1", ... % Stage 01
                       "m_n2"       , "b_2"     , "d_2" , "L_2", ... % Stage 02
                       "m_n3"       , "b_3"     , "d_3" , "L_3", ... % Stage 03
               ... %   Mass         , M. M. Iner, Mass  , M. M. Iner,
                       "m_R"        , "J_R"     , "m_G" , "J_G", ... % Rotor and Generator
               ... %   Length       , Diameter
                       "d_s"        , "L_s"];                        % Main shaft
            
            gamma_full = containers.Map(key_set, ones(18, 1));
            
            switch(aspect)
                case "all"
                    if(numel(gamma) ~= 18)
                        error("gamma must have 18 elements.");
                    end
                    
                    gamma_full = containers.Map(key_set, gamma);
                    
                case "stage"
                    % scaled by stage (i.e. one scale factor for both the
                    % normal module and face width per stage).
                    
                    if(numel(gamma) ~= 3)
                        error("gamma must have 3 elements.");
                    end
                    
                    sub_key ={"m_n1", "b_1", ...
                              "m_n2", "b_2", ...
                              "m_n3", "b_3"};
                          
                    for idx = 1:3
                        gamma_full(sub_key{2*idx - 1}) = gamma(idx, :);
                        gamma_full(sub_key{2*idx})     = gamma(idx, :);
                    end

                case "gear"
                    % scales only the gear dimensions (normal module and 
                    % face width).
                    
                    if(numel(gamma) ~= 6)
                        error("gamma must have 6 elements.");
                    end

                    sub_key ={"m_n1", "b_1", ...
                              "m_n2", "b_2", ...
                              "m_n3", "b_3"};
                          
                    for idx = 1:6
                        gamma_full(sub_key{idx}) = gamma(idx, :);
                    end

                case "K_MMI"
                    % scales only parameters related to the shaft's 
                    % stiffness (length) and mass moment of inertia of 
                    % rotor and generator.

                    if(numel(gamma) ~= 2)
                        error("gamma must have 2 elements.");
                    end
                    
                    sub_key = {"L_1", "L_2", "L_3", ...
                               "J_R", "J_G", ...
                               "L_s"};
                          
                    for idx = 1:3
                        gamma_full(sub_key{idx}) = gamma(1, :);
                    end
                    
                    for idx = 4:5
                        gamma_full(sub_key{idx}) = gamma(2, :);
                    end
                    
                    gamma_full(sub_key{6}) = gamma(1, :);
                    
                case "K_MMI_det"
                    % scales only parameters related to the shaft's 
                    % stiffness (length) and mass moment of inertia of 
                    % rotor and generator.

                    if(numel(gamma) ~= 5)
                        error("gamma must have 6 elements.");
                    end
                    
                    sub_key = {"L_1", "L_2", "L_3", ...
                               "J_R", "J_G", ...
                               "L_s"};
                          
                    for idx = 1:5
                        gamma_full(sub_key{idx}) = gamma(idx, :);
                    end
                    
                    gamma_full(sub_key{end}) = gamma(3, :);
                    
                case "K_MMI_stage"
                    %scales parameters related to the stiffness and mass
                    % moment of inertia, including the gear stages
                    % min |1 - f(sca)/f(ref)|^2 = f(x)
                    % subjected to: S(min) - S(sca) <= 0 = g(x)
                    %               S(ref) - S(sca)  = 0 = h(x)
                    
%                     sub_key = {"m_n1", "b_1", ...
%                                "m_n2", "b_2", ...
%                                "m_n3", "b_3", ...
%                                "L_1" , "L_2", "L_3", ...
%                                "J_R", "J_G", ...
%                                "L_s"};
%                     for idx = 1:3
%                         gamma_full()
%                     end
                    
%                 case "K_MMI_gear"
%                 case "M_MMI_stage_det"
%                 case "K_MMI_gear_det"
                case "K_shaft"
                    % scales only parameters related to the shaft's 
                    % stiffness (length).
                    
                    if(numel(gamma) ~= 3)
                        error("gamma must have 3 elements.");
                    end
                    
                    idx_L = find(contains(key_set, "L"));
                    
                    for idx = 1:3
                        gamma_full(key_set{idx_L(idx)}) = gamma(idx, :);
                    end

                case "K"
                    % scales only stiffness-related parameters. The mesh
                    % stiffness is assumed to be proportinal to the gear's 
                    % face width.
                    
                    if(numel(gamma) ~= 6)
                        error("gamma must have 6 elements.");
                    end

                    gamma_full(2,  :) = gamma(1, :); % face width, stage 01 
                    gamma_full(3,  :) = gamma(2, :); % shaft length
                    gamma_full(6,  :) = gamma(3, :); % face width, stage 02
                    gamma_full(7,  :) = gamma(4, :); % shaft length
                    gamma_full(10, :) = gamma(5, :); % face width, stage 03
                    gamma_full(11, :) = gamma(6, :); % shaft length

                otherwise
                    error("Option [%s] is NOT valid.", aspect);
            end
            
            obj_sca = obj_ref.scale_all(gamma_P, gamma_n, gamma_full);
            
        end
        
        function [obj_sca, gamma_val, res, gamma_sep] = scaled_version(obj_ref, P_scale, n_R_scale, normalize_freq, N_freq, aspect_set, varargin)
            %SCALED_VERSION returns an scaled NREL_5MW object together
            % with the scaling factor for its parameters and residual error
            % from the scaling optimization process. Different scaling can
            % be obtained by using:
            % - model_approach: Different model approaches;
            % - normalize_freq: Normalizing or not the target resonances;
            % - N_freq: choosing how many resonances should be considered
            % during the scaling optimization process;
            % Additionaly, one can define an initial value for gamma to be
            % used on the optimization process
            
            if(~isempty(varargin))
                gamma_prev = varargin{1};
            else
                gamma_prev = ones(18, 1)*0.5;
            end
            
            % Input scaling factors:
            gamma_P = P_scale/obj_ref.P_rated;
            gamma_n = n_R_scale/obj_ref.n_rotor;
            
            % Scaling factors from dimensional analysis:
            gamma_length = nthroot(gamma_P,     3.0);
            gamma_MMI    = power(  gamma_P, 5.0/3.0);
            
            S_H_ref = obj_ref.S_H;
            
            id_1 = "prog:input";
            id_2 = "MATLAB:nearlySingularMatrix";
            warning("off", id_1);
            warning("off", id_2);
            
            fprintf("Scaling NREL 5MW drivetrain with rated power %.1f kW, which is %.2f %% of its reference.\n", gamma_P*[obj_ref.P_rated 100.0]);
            
            %% 1. Stage scaling:
            aspect_1 = "stage";
            
            if(any(aspect_set == aspect_1))
                fprintf("Optimizing drivetrain w.r.t. [%s]...\n", upper(aspect_1));
                
                gamma_0 = mean(gamma_prev([1 2 5 6 9 10]));
                gamma_1 = zeros(obj_ref.N_stg, 1);
                res_1   = zeros(obj_ref.N_stg, 1);

                for idx = 1:obj_ref.N_stg
                    jdx = 2*idx + (-1:0);
                    n_stage = obj_ref.n_out(idx)*gamma_n;

                    [~, gamma_1(idx, :), res_1(idx, :)] = obj_ref.stage(idx).scaled_version(P_scale, n_stage, S_H_ref(jdx), aspect_1, gamma_0);
                    gamma_0 = gamma_1(idx, :);
                end
            elseif(any(aspect_set == "DA"))
                m_n_tmp =[obj_ref.stage.m_n]';
                gamma_1 = Rack.module(m_n_tmp*gamma_length, "calc", "nearest")./m_n_tmp;
                
                res_1 = Inf(3, 1);
            else
                gamma_1 = [mean(gamma_prev([1  2]));
                           mean(gamma_prev([5  6]));
                           mean(gamma_prev([9 10]))];
                       
                res_1 = Inf(3, 1);
            end
            
            gamma_1 = gamma_1([1 1 2 2 3 3]);
            res_1   = res_1  ([1 1 2 2 3 3]);
            
            %% 2. Gears:
            aspect_2 = "gear";
            
            if(any(aspect_set == aspect_2))
                fprintf("Optimizing drivetrain w.r.t. [%s]...\n", upper(aspect_2));
                
                gamma_0 = gamma_1([1 1]);
                gamma_2 = zeros(2*obj_ref.N_stg, 1);
                res_2   = zeros(2*obj_ref.N_stg, 1);

                for idx = 1:obj_ref.N_stg
                    jdx = 2*idx + (-1:0);
                    n_stage = obj_ref.n_out(idx)*gamma_n;

                    [~, gamma_2(jdx, :), res_2(jdx, :)] = obj_ref.stage(idx).scaled_version(P_scale, n_stage, S_H_ref(jdx), aspect_2, gamma_0);
                    gamma_0 = gamma_2(jdx, :);
                end
            elseif(any(aspect_set == "DA"))
                m_n_tmp = [obj_ref.stage.m_n]';
                gamma_mn = Rack.module(m_n_tmp*gamma_length, "calc", "nearest")./m_n_tmp;
                
                gamma_2 = ones(6, 1)*gamma_length;
                gamma_2(1:2:end) = gamma_mn;
                
                res_2 = Inf(6, 1);
            else
                gamma_2 = gamma_prev([1 2 5 6 9 10]);
                res_2 = Inf(6, 1);
            end
            
            if(all(isinf([res_1' res_2'])))
                gamma_12 = gamma_1;
                res_12 = res_1;
            else
                idx_min = res_1 <= res_2;

                gamma_12 = diag(idx_min)*gamma_1 + diag(~idx_min)*gamma_2;
                res_12   = diag(idx_min)*res_1   + diag(~idx_min)*res_2;
            end
            
            obj_12 = obj_ref.scale_aspect(gamma_P, gamma_n, gamma_12, aspect_2);
            
            %% 3. Shaft stiffness and Mass moment of inertia of rotor and generator:
            aspect_3 = "K_MMI";
            
            opt_solver = optimoptions("fmincon", "display", "notify");

            if(any(aspect_set == aspect_3))
                fprintf("Optimizing drivetrain w.r.t. [%s]...\n", upper(aspect_3));
                
                f_n_ref = obj_ref.resonances(N_freq, normalize_freq);

                % function for aspect_3: gamma_P and gamma_n are set to 1.0
                % because the drivetrain is already scaled for these
                % parameters.
                fun_asp  = @(x)(1.0 - obj_12.scale_resonances(N_freq, normalize_freq, 1.0, 1.0, x, aspect_3)./f_n_ref);
                fun_min = @(x)(norm(fun_asp(x))^2);

                gamma_min = ones(2, 1)*1.0e-6;
                gamma_Max = ones(2, 1);

                gamma_0 = [mean(gamma_prev([4 8 12 18])); ... % length
                           mean(gamma_prev([14 16]))]; % mass mom. inertia

                constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities

                [gamma_3, res_3, ~] = fmincon(fun_min, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            elseif(any(aspect_set == "DA"))
                gamma_3 = [gamma_length;
                           gamma_MMI];
                res_3 = Inf;
            else
                gamma_3 = [mean(gamma_prev([4 8 12 18]));
                           mean(gamma_prev([    14 16]))];
                       
                res_3 = Inf;
            end
            
            gamma_3 = gamma_3([1 1 1 2 2]);

            %% 4. Detailed version of 3:
            aspect_4 = "K_MMI_det";
            
            if(any(aspect_set == aspect_4))
                fprintf("Optimizing drivetrain w.r.t. [%s]...\n", upper(aspect_4));
                
                f_n_ref = obj_ref.resonances(N_freq, normalize_freq);

                % function for aspect_4: gamma_P and gamma_n are set to 1.0
                % because the drivetrain is already scaled for these
                % parameters.
                fun_asp  = @(x)(1.0 - obj_12.scale_resonances(N_freq, normalize_freq, 1.0, 1.0, x, aspect_4)./f_n_ref);
                fun_min = @(x)(norm(fun_asp(x))^2);

                gamma_min = ones(5, 1)*1.0e-6;
                gamma_Max = ones(5, 1);

                gamma_0 = gamma_3;

                constraint_fun = @(x)deal([], fun_asp(x)); % inequalities, equalities

                [gamma_4, res_4, ~] = fmincon(fun_min, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
            elseif(any(aspect_set == "DA"))
                gamma_4      = gamma_length*ones(5, 1);
                gamma_4(4:5) = gamma_MMI; 
                
                res_4 = Inf;
            else
                gamma_4 = gamma_prev([4 8 12 14 16]);
                
                res_4 = Inf;
            end
            
            if(all(isinf([res_3 res_4])))
                gamma_34 = gamma_3;
                res_34 = res_3;
            elseif(res_3 <= res_4)
                gamma_34 = gamma_3;
                res_34 = res_3;
            else
                gamma_34 = gamma_4;
                res_34 = res_4;
            end
            
            %% 5. Drivetrain scaling:
            aspect_5 = "Drivetrain";
            
            if(any(aspect_set == aspect_5))
                fprintf("Optimizing drivetrain w.r.t. [%s]...\n", upper(aspect_5));
                
                fprintf("\n\tto be done later...\n\n");
            else
                
            end
            
            %% Post processing:
            % Analysis of the residuals at each step:
            res_pp = [mean(res_1), mean(res_2),    res_3,    res_4];
            idx = ~isinf(res_pp);
            res_pp = res_pp(idx);
            aspect_set = aspect_set(idx);
            
            fprintf("Scale: %.1f kW = %.2f %% of Ref.\n", gamma_P*obj_ref.P_rated, ...
                                                          gamma_P*100.0);

            if(~isempty(res_pp))
                [~, sorted_idx] = sort(res_pp);
                
                fprintf("\tOptimization residua:\n");

                for idx = sorted_idx
                    fprintf("\t%d. %s:\t%.5e\n", idx, aspect_set(idx), res_pp(idx));
                end
            else
                fprintf("\tScaled version obtained by dimensional analysis scaling rules.\n");
            end
            
            gamma_d = nthroot(gamma_P/gamma_n, 3.0);
            
            gamma_val = [gamma_12(1), gamma_12(2), gamma_d   , gamma_34(1), ... % Stage 01
                         gamma_12(3), gamma_12(4), gamma_d   , gamma_34(2), ... % Stage 02
                         gamma_12(5), gamma_12(6), gamma_d   , gamma_34(3), ... % Stage 03
               ... %     Mass       , M. M. Iner , Mass      , M. M. Iner,
                         1.0        , gamma_34(4), 1.0       , gamma_34(5), ... % Rotor and Generator
               ... %     Diameter   , Length
                         gamma_d    , gamma_34(3)]';                            % Main shaft
            
            key_set = ["m_n1"       , "b_1"     , "d_1"   , "L_1" , ... % Stage 01
                       "m_n2"       , "b_2"     , "d_2"   , "L_2" , ... % Stage 02
                       "m_n3"       , "b_3"     , "d_3"   , "L_3" , ... % Stage 03
               ... %   Mass         , M. M. Iner, Mass    , M. M. Iner,
                       "m_R"        , "J_R"     , "m_G" , "J_G", ...    % Rotor and Generator
               ... %   Diameter     , Length
                       "d_s"        , "L_s"];                           % Main shaft
            
            gamma = containers.Map(key_set, gamma_val);
            
            res.stage = res_1;
            res.gear  = res_2;
            res.KJ    = res_3;
            res.KJ2   = res_4;
            res.SG    = res_12;
            res.KJg   = res_34;
            
            gamma_sep.stage = gamma_1;
            gamma_sep.gear  = gamma_2;
            gamma_sep.KJ    = gamma_3;
            gamma_sep.KJ2   = gamma_4;
            
            obj_sca = obj_ref.scale_all(gamma_P, gamma_n, gamma);
            
            warning("on", id_1);
            warning("on", id_2);
            
        end
        
        function [gamma, res, SH, f_n, mode_shape, k_mesh, gamma_asp] = scaled_sweep(obj_ref, P_scale, n_R_scale, normalize_freq, N_freq, aspect_set)
            %SCALED_SWEEP performs a sweep on the rated power parameter of
            % the NREL_5MW object. Returns the scaling factors gamma
            %
            % see also SCALED_VERSION
            %
            
            fprintf("Reference NREL 5MW drivetrain with rated power %.1f kW.\n", obj_ref.P_rated);
            
            n_P = numel(P_scale);
            n_fn = numel(obj_ref.modal_analysis);
            
            gamma = zeros(18, n_P);
            res = struct;
            gamma_asp = struct;
            SH = zeros(numel(obj_ref.S_H), n_P);
            f_n = zeros(n_fn, n_P);
            mode_shape = zeros(n_fn, n_fn, n_P);
            k_mesh = zeros(numel(obj_ref.stage), n_P);
            
            gamma_P = P_scale./obj_ref.P_rated;
            
            [~, MS_ref] = obj_ref.modal_analysis;
            
            SH_ref = obj_ref.S_H;
            SS = [SH_ref, SH_ref];
            
            plot_prop1 = {'lineStyle', '-' , 'lineWidth', 2.0, 'color', [346.6667e-3   536.0000e-3   690.6667e-3]};
            plot_prop2 = {'lineStyle', '--', 'lineWidth', 2.0, 'color', [915.2941e-3   281.5686e-3   287.8431e-3]};
            
            figure("units", "centimeters", "position", [5.0 5.0 34.0 12.0]);
            subplot(2, 6, 1:3)
            rectangle(obj_ref);
            xlim([0 6000])
            ylim([-1 1]*1500)
            title(sprintf("Reference: %.1f kW", obj_ref.P_rated));
            
            subplot(2, 6, 7:9)
            axis equal;
            xlim([0 6000])
            ylim([-1 1]*1500)
            
            for idx = 1:3
                subplot(2, 6, idx + 3)
                plot(1:14, MS_ref(:, idx)*90.0, plot_prop1{:});
                
                if(idx ~= 3)
                    subplot(2, 6, idx + 9)
                    plot(1:14, MS_ref(:, idx + 3)*90.0, plot_prop1{:})
                else
                    subplot(2, 6, 12)
%                     bar(SS);
                end
            end
            
            tick_x = ["R"  , ...
                      "c_1", "p_{11}", "p_{12}", "p_{13}", "s_1", ...
                      "c_2", "p_{21}", "p_{22}", "p_{23}", "s_2", ...
                      "W_3", "P_3", ...
                      "G"];
            idx_tick = 2:2:14;
            
            font_setting  = {"fontName", "Times", "fontSize", 12.0};
            latex_setting = {"fontName", "Times", "fontSize", 12.0, "interpreter", "LaTeX"};
            
            fig_axes = findobj(gcf, "Type", "Axes");
            
            set(get(fig_axes(1), "ylabel"), "string", "$S_H$"             , latex_setting{:});
            set(get(fig_axes(5), "ylabel"), "string", "$\theta_i$, [deg.]", latex_setting{:});
            set(get(fig_axes(6), "ylabel"), "string", "$\theta_i$, [deg.]", latex_setting{:});

            set(fig_axes(2:6), 'xlim'      , [1 14]);
            set(fig_axes(2:6), 'ylim'      , [-1 1]*100);
            set(fig_axes(2:6), 'xtick'     ,        idx_tick);
            set(fig_axes(2:6), 'xticklabel', tick_x(idx_tick));
            set(fig_axes     , font_setting{:});
            
            label_figure(fig_axes, font_setting);
            
            fig_name = @(i)(sprintf("plots\\sweep_scale\\scale_%d_%d", i, n_P));
            
            file_name = "plots\sweep_scale\scale_sweep.gif";
            
            k_P = mean(diff(P_scale));
            
            if(k_P > 0.0)
                gamma_0 = ones(18, 1)*1.0e-3;
            else
                gamma_0 = ones(18, 1);
            end
            

            for idx = 1:n_P
                [obj_sca, gamma_0, res_idx, gamma_idx] = obj_ref.scaled_version(P_scale(idx), n_R_scale, normalize_freq, N_freq, aspect_set, gamma_0);
                
                gamma(:, idx) = gamma_0;
                res(idx).stage = res_idx.stage;
                res(idx).gear  = res_idx.gear;
                res(idx).KJ    = res_idx.KJ;
                res(idx).KJ2   = res_idx.KJ2;
                
                gamma_asp(idx).stage = gamma_idx.stage;
                gamma_asp(idx).gear  = gamma_idx.gear;
                gamma_asp(idx).KJ    = gamma_idx.KJ;
                gamma_asp(idx).KJ2   = gamma_idx.KJ2;
                
                SH(:, idx) = obj_sca.S_H;
                k_mesh(:, idx) = [obj_sca.stage.k_mesh]';
                [f_n(:, idx), mode_shape(:, :, idx)] = obj_sca.modal_analysis;

                subplot(2, 6, 7:9);
                cla;
                rectangle(obj_sca);
                xlim([0 6000]);
                ylim([-1 1]*1500);
                title(sprintf("Scale: %.1f kW = %.2f %% of Ref.", obj_sca.P_rated, gamma_P(idx)*100.0));
                SS(:, 2) = SH(:, idx);
                
                for jdx = 1:3
                    subplot(2, 6, jdx + 3);
                    cla;
                    hold on;
                    plot(1:14, MS_ref(:, jdx)*90.0, plot_prop1{:});
                    plot(1:14, mode_shape(:, jdx, idx)*90.0, plot_prop2{:});

                    if(jdx ~= 3)
                        subplot(2, 6, jdx + 9);
                        cla;
                        hold on;
                        plot(1:14, MS_ref(:, jdx + 3)*90.0, plot_prop1{:})
                        plot(1:14, mode_shape(:, jdx + 3, idx)*90.0, plot_prop2{:});
                    else
                        subplot(2, 6, 12);
                        cla;
                        b = bar(SS);
                        yline(1.25, "color", "k");
                        ylim([1.0 2.0]);
                        yticks(1.0:0.25:2.0);
                        b(1).FaceColor = plot_prop1{end};
                        b(2).FaceColor = plot_prop2{end};
                    end
                end

                set(get(fig_axes(1), "ylabel"), "string", "$S_H$, [-]"        , latex_setting{:});
                set(fig_axes(1)  , 'xticklabel', ["s_1", "p_{1i}", "s_2", "p_{2i}", "P_3", "W_3"]);
                
                set(fig_axes, font_setting{:});
                savefig(gcf, fig_name(idx));
                print(fig_name(idx), '-dpng');
                saveasGIF(file_name, idx);
                
            end
            
        end
        
        
    end
    
end
