classdef NREL_5MW < Drivetrain
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
    % see also DRIVETRAIN, DTU_10MW.
    
    properties
        gamma_P (1, 1)                {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % Scaling factor for the rated power
        gamma_n (1, 1)                {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0;    % Scaling factor for the rotor speed
    end
    
    methods
        function obj = NREL_5MW(varargin)
            gamma = {'P'   , 1.0, ...
                     'n'   , 1.0, ...
                     'm_n1', 1.0, ...
                     'b_1' , 1.0, ...
                     'd_1' , 1.0, ...
                     'L_1' , 1.0, ...
                     'm_n2', 1.0, ...
                     'b_2' , 1.0, ...
                     'd_2' , 1.0, ...
                     'L_2' , 1.0, ...
                     'm_n3', 1.0, ...
                     'b_3' , 1.0, ...
                     'd_3' , 1.0, ...
                     'L_3' , 1.0, ...
                     'J_R' , 1.0, ...
                     'J_G' , 1.0, ...
                     'd_S' , 1.0, ...
                     'L_S' , 1.0};
            
            gamma = process_varargin(gamma, varargin{:});
            gamma = scaling_factor(gamma);
                        
            N_st = 3;
            stage = [Gear_Set Gear_Set Gear_Set];
            
            for idx = 1:N_st
                stg = NREL_5MW.gear_stage(idx);
                
                gamma_idx = gamma.ends_with(num2str(idx));
                stage(idx) = stg.scale_by(gamma_idx);
            end
            
            P_r = gamma('P')* 5.0e3; % [kW], Rated power
            n_r = gamma('n')*12.1;   % [1/min.], Input speed
            
            LSS = NREL_5MW.shaft(0);
            
            inp_shaft = Shaft(LSS.d*gamma("d_S"), ...
                              LSS.L*gamma("L_S"));
            
            m_R = 110.0e3;
            J_R = 57231535.0  *gamma("J_R"); % according to [1, 3]
            m_G =     1900.0;
            J_G =      534.116*gamma("J_G");
            
            obj@Drivetrain('N_stage',    N_st, ...
                           'stage',      stage, ...
                           'P_rated',    P_r, ...
                           'n_rotor',    n_r, ...
                           'main_shaft', inp_shaft, ...
                           'm_Rotor',    m_R, ...
                           'J_Rotor',    J_R, ...
                           'm_Gen',      m_G, ...
                           'J_Gen',      J_G, ...
                           'S_Hmin',     1.25, ...
                           'S_Fmin',     1.56, ...
... %                            'dynamic_model', @Kahraman_94);
                           'dynamic_model', @Lin_Parker_99);
                           
            [obj.S_H_val, obj.S_F_val, obj.S_shaft_val] = obj.safety_factors();

            obj.gamma = scaling_factor(gamma.name, ones(length(gamma), 1));
        end
    end
    
    %% Calculation:
    methods
        function [SH, SF, SShaft] = safety_factor_stage(obj, idx)
            L_h    = 20*365*24;
            K_A    = 1.25;
            
            calc = ISO_6336(obj.stage(idx), 'calculation', 'KISSsoft', ...
                                            'P_rated'    , obj.P_rated, ...
                                            'n_out'      , obj.n_out(idx), ...
                                            'S_Hmin'     , obj.S_Hmin, ...
                                            'S_Fmin'     , obj.S_Fmin, ...
                                            'L_h'        , L_h, ... % [h], Required life
                                            'K_A'        , K_A);    % [-], Application factor
            
            [SH, SF] = calc.safety_factors('lubricant_ID', '11220', ...
                                           'nu_40'       , 220.0, ...
                                           'stage_idx'   ,   0, ...
                                           'save_report' , false, ...
                                           'show_report' , false, ...
                                           'line'        ,   4);

            S_ut = Material.S_ut*1.0e-6; % [MPa], Tensile strength
            S_y  = Material.S_y*1.0e-6;  % [MPa], Yield strength
            K_f  = 1.0;                  % [-],   Fatigue stress-concentration factor for bending
            K_fs = 1.0;                  % [-],   Fatigue stress-concentration factor for torsion
            
            SShaft = obj.stage(idx).output_shaft.safety_factors(S_ut, S_y, K_f, K_fs, obj.T_out(idx));
        end
        
    end
    
    methods(Static)
        function brg = bearing(idx)
            %BEARING returns the bearing sets of each stage of the NREL 
            % 5 MW wind turbine drivetrain according to [1], tables VI and
            % IX.
            %
            cx = 453.0;         cy = 42.0e3;            cz = 30600.0;
                                cb = 34.3;              cg = 47.8;
            switch(idx)
                case 0 % Main shaft
                    INP_A     = Bearing('name', 'INP_A', 'type'  ,  'CARB', ...
                                        'K_x' ,   0.0e4, 'K_y'   ,  1.5e10, 'K_z'    , 1.5e10, ...
                                                         'K_beta',   5.0e6, 'K_gamma',  5.0e6, ...
                                        'C_x' ,  0.0*cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'OD'  ,  1750.0, 'ID'    ,  1250.0, 'B'      ,  375.0);
                    INP_B     = Bearing('name', 'INP_B', 'type'  , 'SRB'  , ...
                                        'K_x' ,  4.06e8, 'K_y'   , 1.54e10, 'K_z'    , 1.54e10, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'OD'  ,  1220.0, 'ID'    ,   750.0, 'B'      ,   365.0);
                    
                    brg = [INP_A, INP_B];
                    
                case 1
                    PL_A      = Bearing('name',  'PL_A', 'type'  ,  'CRB', ...
                                        'K_x' ,   9.1e4, 'K_y'   ,  9.4e9, 'K_z'    , 3.2e9, ...
                                                         'K_beta',  1.4e6, 'K_gamma', 4.5e6, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'OD'  ,   600.0, 'ID'    ,  400.0, 'B'      , 272.0);
                    PL_B      = Bearing('name',  'PL_B', 'type'  ,  'CRB', ...
                                        'K_x',    9.1e4, 'K_y'   ,  9.4e9, 'K_z'    , 3.2e9, ...
                                                         'K_beta',  1.4e6, 'K_gamma', 4.5e6, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'OD' ,    600.0, 'ID'    ,  400.0, 'B'      , 272.0);
                    PLC_A     = Bearing('name', 'PLC_A', 'type'  ,  'SRB', ...
                                        'K_x',    6.6e4, 'K_y'   ,  1.7e9, 'K_z'    , 1.1e9, ...
                                                         'K_beta',  5.6e5, 'K_gamma', 1.3e5, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'OD' ,   1030.0, 'ID'    ,  710.0, 'B'      , 315.0);
                    PLC_B     = Bearing('name', 'PLC_B', 'type'  ,  'CRB', ...
                                        'K_x',    6.6e7, 'K_y'   ,  1.7e9, 'K_z'    , 1.1e9, ...
                                                         'K_beta',  5.6e5, 'K_gamma', 1.3e5, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'OD' ,   1220.0, 'ID'    , 1000.0, 'B'      , 128.0);
                    
                    brg = [PL_A,  PL_B,            ... % Planet
                           PLC_A, PLC_B,           ... % Carrier
                           Bearing('name', 'sun'), ... % Sun
                           Bearing('name', 'ring')];   % Ring
                           
                case 2
                    IMS_PL_A  = Bearing('name',  "IMS_PL_A", 'type'  , "CRB" , ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'OD'  ,       520.0, 'ID'    ,  380.0, 'B'      , 140.0);
                    IMS_PL_B  = Bearing('name',  "IMS_PL_B", 'type'  , "CRB" , ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'OD'  ,       520.0, 'ID'    ,  380.0, 'B'      , 140.0);
                    IMS_PLC_A = Bearing('name', "IMS_PLC_A", 'type'  , "CARB", ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'OD'  ,      1030.0, 'ID'    ,  710.0, 'B'      , 236.0);
                    IMS_PLC_B = Bearing('name', "IMS_PLC_B", 'type'  , "CRB" , ...
                                        'K_x' ,       9.1e7, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'OD'  ,       870.0, 'ID'    ,  600.0, 'B'      , 155.0);
                    
                    brg = [IMS_PL_A,  IMS_PL_B,    ... % Planet
                           IMS_PLC_A, IMS_PLC_B,   ... % Carrier
                           Bearing('name', 'sun'), ... % Sun
                           Bearing('name', 'ring')];   % Ring
                    
                case 3
                    IMS_A     = Bearing('name', "IMS_A", 'type'  , "CRB", ...
                                        'K_x' ,     0.0, 'K_y'   , 6.0e7, 'K_z'    , 1.2e9, ...
                                                         'K_beta', 7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   360.0, 'ID'    , 200.0, 'B'      ,  98.0);
                    IMS_B     = Bearing('name', "IMS_B", 'type'  , "TRB", ...
                                        'K_x' ,   7.4e7, 'K_y'   , 5.0e8, 'K_z'    , 5.0e8, ...
                                                         'K_beta', 1.6e6, 'K_gamma', 1.8e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   460.0, 'ID'    , 200.0, 'B'      , 100.0);
                    IMS_C     = Bearing('name', "IMS_C", 'type'  , "TRB", ...
                                        'K_x' ,   7.8e7, 'K_y'   , 7.4e8, 'K_z'    , 3.3e8, ...
                                                         'K_beta', 1.1e6, 'K_gamma', 2.5e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   460.0, 'ID'    , 200.0, 'B'      , 100.0);
                    HS_A      = Bearing('name',  "HS_A", 'type'  , "CRB", ...
                                        'K_x' ,   1.3e8, 'K_y'   , 8.2e8, 'K_z'    , 8.2e8, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   500.0, 'ID'    , 400.0, 'B'      , 100.0);
                    HS_B      = Bearing('name',  "HS_B", 'type'  , "TRB", ...
                                        'K_x' ,   6.7e7, 'K_y'   , 8.0e8, 'K_z'    , 1.3e8, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   550.0, 'ID'    , 410.0, 'B'      ,  86.0);
                    HS_C      = Bearing('name',  "HS_C", 'type'  , "TRB", ...
                                        'K_x' ,   8.0e7, 'K_y'   , 1.0e9, 'K_z'    , 7.3e7, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   550.0, 'ID'    , 410.0, 'B'      ,  86.0);
                    
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
                    d = 700.0;
                    L = 2000.0;

                case 1 % ISS
%                     d_1 = 800.0;
%                     L_1 = 750.0;
                    d = 533.0;
                    L = 500.0;

                case 2 % HS-IS
%                     d_2 = 500.0;
%                     L_2 = 1000.0;
                    d = 333.0;
                    L = 666.0;

                case 3 % HSS
%                     d_3 = 500.0;
%                     L_3 = 1500.0;
                    d = 333.0;
                    L = 1000.0;

                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", idx);
            end
            
            sft = Shaft(d, L);
        end
        
        function g_set = gear_stage(idx)
            %GEAR_STAGE returns the gear stages of the NREL 5 MW wind
            % turbine drivetrain according to [4], Table V. The values for
            % the tip alteration coefficients were taken from KISSsoft.
            %

            Q = 6.0;
            Ra = 0.8;
            
            switch(idx)
                case 1
                    p    =   3;         % [-],    Number of planet gears
                    m_n  =  45.0;       % [mm],   Normal module
                    beta =   0.0;       % [deg.], Helix angle (at reference cylinder)
                    b    = 491.0;       % [mm],   Face width
                    a_w  = 863.0;       % [mm],   Center distance
                    z_s  =  19;         % [-],    Number of teeth (sun)    [WHEEL]
                    z_p  =  17;         % [-],    Number of teeth (planet) [PINION]
                    z_r  =  56;         % [-],    Number of teeth (ring)
                    x_s  =   0.617;     % [-],    Profile shift coefficient (sun)
                    x_p  =   0.802;     % [-],    Profile shift coefficient (planet)
                    x_r  =  -0.501;     % [-],    Profile shift coefficient (ring)
                    k_s  = -10.861/m_n; % [-],    Tip alteration coefficient (sun)
                    k_p  = -10.861/m_n; % [-],    Tip alteration coefficient (planet)
                    k_r  =   0.0;       % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs = 80.0/171.0;
                    bore_Rp = 80.0/153.0;
                    bore_Rr =  1.2;
                    
                    z = [z_s z_p z_r];
                    x = [x_s x_p x_r];
                    k = [k_s k_p k_r];
                    bore_ratio = [bore_Rs bore_Rp bore_Rr];
                    config = "planetary";
                    
                case 2
                    p    =   3;        % [-],    Number of planet gears
                    m_n  =  21.0;      % [mm],   Normal module
                    beta =   0.0;      % [deg.], Helix angle (at reference cylinder)
                    b    = 550.0;      % [mm],   Face width
                    a_w  = 584.0;      % [mm],   Center distance
                    z_s  =  18;        % [-],    Number of teeth (sun)    [PINION]
                    z_p  =  36;        % [-],    Number of teeth (planet) [WHEEL]
                    z_r  =  93;        % [-],    Number of teeth (ring)
                    x_s  =   0.389;    % [-],    Profile shift coefficient (sun)
                    x_p  =   0.504;    % [-],    Profile shift coefficient (planet)
                    x_r  =   0.117;    % [-],    Profile shift coefficient (ring)
                    k_s  =  -1.75/m_n; % [-],    Tip alteration coefficient (sun)
                    k_p  =  -1.75/m_n; % [-],    Tip alteration coefficient (planet)
                    k_r  =   0.0;      % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs = 100.0/189.0;
                    bore_Rp =  95.0/189.0;
                    bore_Rr =   1.2;

                    z = [z_s z_p z_r];
                    x = [x_s x_p x_r];
                    k = [k_s k_p k_r];
                    bore_ratio = [bore_Rs bore_Rp bore_Rr];
                    config = "planetary";

                case 3
                    p    =   1;         % [-],    Number of planet gears
                    m_n  =  14.0;       % [mm],   Normal module
                    beta =  10.0;       % [deg.], Helix angle (at reference cylinder)
                    b    = 360.0;       % [mm],   Face width
                    a_w  = 861.0;       % [mm],   Center distance
                    z_1  =  24;         % [-],    Number of teeth (pinion)
                    z_2  =  95;         % [-],    Number of teeth (wheel)
                    x_1  =   0.480;     % [-],    Profile shift coefficient (pinion)
                    x_2  =   0.669;     % [-],    Profile shift coefficient (wheel)
                    k_1  =  -0.938/m_n; % [-],    Tip alteration coefficient (pinion)
                    k_2  =  -0.938/m_n; % [-],    Tip alteration coefficient (wheel)
                    
                    bore_R1 = 1809.0/3086.0;
                    bore_R2 = 3385.0/9143.0;

                    z = [z_1 z_2];
                    x = [x_1 x_2];
                    k = [k_1 k_2];
                    bore_ratio = [bore_R1 bore_R2];
                    config = "parallel";
                    
                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", idx);
            end
            
            g_set = Gear_Set('configuration', config, ...
                             'm_n'          , m_n, ...
                             'z'            , z, ...
                             'b'            , b, ...
                             'x'            , x, ...
                             'beta'         , beta, ...
                             'k'            , k, ...
                             'bore_ratio'   , bore_ratio, ...
                             'N_p'          , p, ...
                             'a_w'          , a_w, ...
                             'bearing'      , NREL_5MW.bearing(idx), ...
                             'shaft'        , NREL_5MW.shaft(idx));
        end
        
        function property_estimation()
            u = 97.0;
            m_r = 110.0e3; % [kg], according to [2], Table 1-1
            J_g = 534.116; % [kg-m^2], according to [2], Table 5-1
            k_eq = 867637.0e3; % [N-m/rad], taken from [2], Table 5-1
            
            freq_free = 2.18; % [Hz], taken from [3], p. 45
            freq_fix  = mean([0.6205 0.6094]); % [Hz], taken from [2], Table 9-1. Mean of the results obtained using FAST and ADAMS, respectively.
            
            J_r_fix   = k_eq/(2*pi*freq_fix)^2;
            J_r_free  = (k_eq*J_g*u^2)/(J_g*(2.0*pi*u*freq_free)^2 - k_eq);
            J_r_ratio = ((freq_free/freq_fix)^2 - 1.0)*J_g*u^2;
            
            R_cyl = sqrt(2.0*J_r_ratio/m_r); % [m], cylinder radius
            h_cyl = m_r/(Material.rho*pi*R_cyl^2); % [m], cylinder height
            
            fprintf("Rotor mass moment of inertia using:\n")
            fprintf("\t Rigid free-fixed resonance: %3.4e [Hz]\t %3.4e [kg-m^2]\n", freq_fix, J_r_fix);
            fprintf("\t Rigid free-free  resonance: %3.4e [Hz]\t %3.4e [kg-m^2]\n", freq_free, J_r_free);
            fprintf("\t Ratio of the frequencies above: \t\t %3.4e [kg-m^2]\n", J_r_ratio);
            fprintf("Cylindric rotor dimensions: R = %.3f \t h = %.3f [m]\n", R_cyl, h_cyl);
            
        end
    end
    
    methods
        update_subvar(obj)
    end
    
end
