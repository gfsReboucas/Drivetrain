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
            
            N_st = 3;
            stage = [Gear_Set Gear_Set Gear_Set];

            %      Gear stage               , Output shaft
            %      Normal module, Face width, Length, Diameter
            par = ["m_n1"       , "b_1"     , "d_1" , "L_1", ... % Stage 01
                   "m_n2"       , "b_2"     , "d_2" , "L_2", ... % Stage 02
                   "m_n3"       , "b_3"     , "d_3" , "L_3", ... % Stage 03
                   "J_R",                                    ... % M. M. Inertia (rotor)
                   "J_G",                                    ... % M. M. Inertia (generator)
                                              "d_s" , "L_s"]'; % Main shaft

            if(isempty(varargin))

                for idx = 1:N_st
                    stage(idx) =  NREL_5MW.gear_stage(idx);
                end

                P_r = 5.0e3; % [kW], Rated power
                n_r = 12.1; % [1/min.], Input speed
                inp_shaft = NREL_5MW.shaft(0);

                m_R = 110.0e3; % [kg], according to [2], Table 1-1
                J_R = 57231535.0; % [kg-m^2], according to [3], p. 45
%                 J_R = 58.1164e6; % according to property_estimation
                m_G = 0.02*m_R; % 1900.0; % [kg], but according to ???
                J_G = 534.116; % [kg-m^2], according to [2], Table 5-1

                gm_val = ones(size(par));
                
            elseif(length(varargin) == 3)
                if(isa(varargin{3}, "scaling_factor"))
                    gm_P = varargin{1};
                    gm_n = varargin{2};
                    gm   = varargin{3};
                    
                    for idx = 1:N_st
                        stg = NREL_5MW.gear_stage(idx);
                        
                        jdx = 4*idx + (-3:0);
                        gm_stg = gm(jdx);
                        stage(idx) = stg.scale_aspect(gm_stg, "Gear_Set");
                    end
                    
                    P_r = gm_P*5.0e3; % [kW], Rated power
                    n_r = gm_n*12.1; % [1/min.], Input speed
                    
                    LSS = NREL_5MW.shaft(0);
                    
                    inp_shaft = Shaft(LSS.d*gm("d_s"), ...
                                      LSS.L*gm("L_s"));
                                  
                    m_R = 110.0e3;
                    J_R = 57231535.0*gm("J_R"); % according to [1, 3]
                    m_G = 0.02*m_R;
                    J_G = 534.116*gm("J_G");
                    
                    gm_val = gm.value;
                end
            end
            
            obj@Drivetrain(N_st, stage, P_r, n_r, inp_shaft, m_R, J_R, m_G, J_G);
            obj.dynamic_model =  "Kahraman_1994";
            
            obj.gamma = scaling_factor(par, gm_val);
            
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

            alpha_n = 20.0;        % [deg.],   Pressure angle (at reference cylinder)
            rack_type = "A";       % [-],      Type of the basic rack from A to D
            Q = 6.0;
            
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
                    bore_R = [bore_Rs bore_Rp bore_Rr];
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
                    bore_R = [bore_Rs bore_Rp bore_Rr];
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
                    bore_R = [bore_R1 bore_R2];
                    config = "parallel";
                    
                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", idx);
            end
            
            g_set = Gear_Set(config, m_n, alpha_n, z, b, x, beta, k, bore_R, p, a_w, rack_type, NREL_5MW.bearing(idx), NREL_5MW.shaft(idx), Q);
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
