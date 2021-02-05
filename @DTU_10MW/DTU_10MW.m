classdef DTU_10MW < Drivetrain
    %DTU_10MW This class contains some of the properties of the DTU 10MW
    % wind turbine gearbox proposed by Wang et. al. [1]. More information
    % regarding the DTU 10MW wind turbine can be found at [2].
    %
    % [1] S. Wang, A. Nejad and T. Moan, "On design, modelling, and
    % analysis of a 10-MW medium-speed drivetrain for offshore wind
    % turbines", Wind Energy, 2020. doi:10.1002/we.2476
    % [2] C. Bak; F. Zahle; R. Bitsche; T. Kim; A. Yde; L.C. Henriksen;
    % P.B. Andersen; A. Natarajan, M.H. Hansen; “Design and performance of
    % a 10 MW wind turbine”, J. Wind Energy, To be accepted.
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
    % see also DRIVETRAIN, NREL_5MW.
    
    properties
        gamma_P (1, 1)                {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % Scaling factor for the rated power
        gamma_n (1, 1)                {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % Scaling factor for the rotor speed
    end
    
    methods
        function obj = DTU_10MW(varargin)
            gamma = {'P'   , 1.0, 'n'   , 1.0, ...
                     'm_n1', 1.0, 'b_1' , 1.0, 'd_1' , 1.0, 'L_1' , 1.0, ...
                     'm_n2', 1.0, 'b_2' , 1.0, 'd_2' , 1.0, 'L_2' , 1.0, ...
                     'm_n3', 1.0, 'b_3' , 1.0, 'd_3' , 1.0, 'L_3' , 1.0, ...
                     'J_R' , 1.0, 'J_G' , 1.0, ...
                     'd_S' , 1.0, 'L_S' , 1.0};
            
            gamma = scaling_factor.process_varargin(gamma, varargin);
            gamma = scaling_factor(gamma);
                        
            N_st = 3;
            stage = repmat(Gear_Set(), 1, N_st);
            
            for idx = 1:N_st
                stg = DTU_10MW.gear_stage(idx);
                
                gamma_idx = gamma.ends_with(num2str(idx));
                stage(idx) = stg.scale_by(gamma_idx);
            end
            
            P_r = gamma('P')* 5.0e3; % [kW], Rated power
            n_r = gamma('n')*12.1;   % [1/min.], Input speed
            
            LSS = DTU_10MW.shaft(0);
            
            inp_shaft = Shaft(LSS.d*gamma('d_S'), ...
                              LSS.L*gamma('L_S'), ...
                              LSS.bearing, ...
                              LSS.material);

            m_R = 227962.0;                   % [kg], according to [2], Table 2.1
            J_R =    156.7374e6*gamma('J_R'); % [kg-m^2] according to property_estimation
            m_G =      0.02*m_R;              % [kg], according to ???
            J_G =   1500.5     *gamma('J_G'); % [kg-m^2] according to [2], Table 6.3
            
            obj@Drivetrain('N_stage'      , N_st, ...
                           'stage'        , stage, ...
                           'P_rated'      , P_r, ...
                           'n_rotor'      , n_r, ...
                           'main_shaft'   , inp_shaft, ...
                           'm_Rotor'      , m_R, ...
                           'J_Rotor'      , J_R, ...
                           'm_Gen'        , m_G, ...
                           'J_Gen'        , J_G, ...
                           'S_Hmin'       , 1.25, ...
                           'S_Fmin'       , 1.56, ...
                           'dynamic_model', @Kahraman_94);

            obj.gamma = scaling_factor(gamma.name, gamma.value);

            [obj.S_H_val, obj.S_F_val, ...
                obj.S_shaft_val] = obj.safety_factors();
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
                    p_1    = 5;            % [-],    Number of planet gears
                    m_n1   =  30.0;        % [mm],   Normal module
                    beta_1 =   8.0;        % [deg.], Helix angle (at reference cylinder)
                    b_1    = 800.0;        % [mm],   Face width
                    a_w1   = 877.033;      % [mm],   Center distance
                    z_s1   =  26;          % [-],    Number of teeth (sun)    [WHEEL]
                    z_p1   =  31;          % [-],    Number of teeth (planet) [PINION]
                    z_r1   =  89;          % [-],    Number of teeth (ring)
                    x_s1   =   0.2702;     % [-],    Profile shift coefficient (sun)
                    x_p1   =   0.2093;     % [-],    Profile shift coefficient (planet)
                    x_r1   =  -0.1591;     % [-],    Profile shift coefficient (ring)
                    k_s1   = -0.756/m_n1; % [-],    Tip alteration coefficient (sun)
                    k_p1   = -0.756/m_n1; % [-],    Tip alteration coefficient (planet)
                    k_r1   =   0.0;        % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs1 = 80.0/171.0;
                    bore_Rp1 = 80.0/153.0;
                    bore_Rr1 = 1.2;
                    
                    z_1 = [z_s1 z_p1 z_r1];
                    x_1 = [x_s1 x_p1 x_r1];
                    k_1 = [k_s1 k_p1 k_r1];
                    bore_R1 = [bore_Rs1 bore_Rp1 bore_Rr1];
%                                Sun,   Planet,   Ring,   Carrier
                    bearing_1 = DTU_10MW.bearing(1);
                    shaft_1 = DTU_10MW.shaft(1);
                    
                    g_set = Gear_Set("planetary", m_n1, alpha_n, z_1, b_1, x_1, beta_1, k_1, bore_R1, p_1, a_w1, rack_type, bearing_1, shaft_1);

                case 2
                    p_2    =   3;        % [-],    Number of planet gears
                    m_n2   =  20.0;      % [mm],   Normal module
                    beta_2 =   8.0;      % [deg.], Helix angle (at reference cylinder)
                    b_2    = 520.0;      % [mm],   Face width
                    a_w2   = 684.273;    % [mm],   Center distance
                    z_s2   =  26;        % [-],    Number of teeth (sun)    [PINION]
                    z_p2   =  41;        % [-],    Number of teeth (planet) [WHEEL]
                    z_r2   = 109;        % [-],    Number of teeth (ring)
                    x_s2   =   0.2787;   % [-],    Profile shift coefficient (sun)
                    x_p2   =   0.1213;   % [-],    Profile shift coefficient (planet)
                    x_r2   =  -0.0024;   % [-],    Profile shift coefficient (ring)
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
                    bearing_2 = DTU_10MW.bearing(2);
                    shaft_2 = DTU_10MW.shaft(2);
                    
                    g_set = Gear_Set("planetary", m_n2, alpha_n, z_2, b_2, x_2, beta_2, k_2, bore_R2, p_2, a_w2, rack_type, bearing_2, shaft_2);

                case 3
                    p_3    = 1;            % [-],    Number of planet gears
                    m_n3   = 18.0;         % [mm],   Normal module
                    beta_3 = 12.0;         % [deg.], Helix angle (at reference cylinder)
                    b_3    = 500.0;        % [mm],   Face width
                    a_w3   = 825.885;        % [mm],   Center distance
                    z_13   =  28;          % [-],    Number of teeth (pinion)
                    z_23   =  61;          % [-],    Number of teeth (wheel)
                    x_13   =  0.2976;       % [-],    Profile shift coefficient (pinion)
                    x_23   =  0.1024;       % [-],    Profile shift coefficient (wheel)
                    k_13   =  -0.938/m_n3; % [-],    Tip alteration coefficient (pinion)
                    k_23   =  -0.938/m_n3; % [-],    Tip alteration coefficient (wheel)
                    
                    bore_R13 = 1809.0/3086.0;
                    bore_R23 = 3385.0/9143.0;

                    z_3 = [z_13 z_23];
                    x_3 = [x_13 x_23];
                    k_3 = [k_13 k_23];
                    bore_R3 = [bore_R13 bore_R23];
%                                 Pinion, Wheel
                    bearing_3 = DTU_10MW.bearing(3);
                    shaft_3 = DTU_10MW.shaft(3);
                    
                    g_set = Gear_Set("parallel", m_n3, alpha_n, z_3, b_3, x_3, beta_3, k_3, bore_R3, p_3, a_w3, rack_type, bearing_3, shaft_3);
                    
                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", idx);
            end
        end
        
        function property_estimation()
            u = 50.0;
            m_r = 227962.0; % [kg], according to [2], Table 2.1
            J_g = 1500.5;        % [kg-m^2], taken from [2], Table 6.3
            k_eq = 2317025352.0; % [N-m/rad], taken from [2], Table 6.3
            
            freq_free = 4.003; % [Hz], taken from [2], Table 6.2
            freq_fix  = 0.612; % [Hz], taken from [2], Table 6.2
            
            rho = Material().rho*1.0e9;
            
            J_r_fix   = k_eq/(2.0*pi*freq_fix)^2;
            J_r_free  = (k_eq*J_g*u^2)/(J_g*(2.0*pi*u*freq_free)^2 - k_eq);
            J_r_ratio = ((freq_free/freq_fix)^2 - 1.0)*J_g*u^2;
            
            R_cyl = sqrt(2.0*J_r_ratio/m_r); % [m], cylinder radius
            h_cyl = m_r/(rho*pi*R_cyl^2); % [m], cylinder height
            
            f = @(x)([(pi*rho/m_r)            .*x(1, :).*(1.0 - x(3, :).^2).*x(2, :)^2 - 1.0; ...
                      (pi*rho/(2.0*J_r_ratio)).*x(1, :).*(1.0 - x(3, :).^4).*x(2, :)^4 - 1.0]);
            x = fmincon(@(x)(norm(f(x)))^2, rand*ones(3,1), [], [], [], [], zeros(3, 1), [Inf, Inf, 1.0]');
            r_ext = x(2);           r_int = x(3)*r_ext;           h = x(1);

            fprintf("Rotor mass moment of inertia using:\n");
            fprintf("\t Rigid free-fixed resonance: %3.4e [Hz]\t %3.4e [kg-m^2]\n", freq_fix, J_r_fix);
            fprintf("\t Rigid free-free  resonance: %3.4e [Hz]\t %3.4e [kg-m^2]\n", freq_free, J_r_free);
            fprintf("\t Ratio of the frequencies above: \t\t %3.4e [kg-m^2]\n", J_r_ratio);
            fprintf("Rotor dimensions assuming:\n");
            fprintf("\t Cylindric geometry: R = %.3f \t h = %.3f [m]\n", R_cyl, h_cyl);
            fprintf("\t Cylindric tube geometry, opt.: R_ext = %.3f \t R_int = %.3f \t h = %.3f [m]\n", r_ext, r_int, h);
            
        end
    end
    
    %% Set methods:
    methods
        function obj = set.gamma_P(obj, val)
            obj.gamma_P = val;
            obj.P_rated = obj.P_rated*obj.gamma_P;
        end
        
        function obj = set.gamma_n(obj, val)
            obj.gamma_n = val;
            obj.n_rotor = obj.n_rotor*obj.gamma_n;
        end
    end
    %% Get methods:
    methods
%         function val = get.gamma(obj)
%             val = obj.gamma;
%         end
    end
    
end
