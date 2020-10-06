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
        gamma_P (1, 1) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0; % Scaling factor for the rated power
        gamma_n (1, 1) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0; % Scaling factor for the rotor speed
    end
    
    methods
        function obj = NREL_5MW(varargin)
            gamma = {'P'   , 1.0, 'n'   , 1.0, ...
                     'm_n1', 1.0, 'b_1' , 1.0, 'd_1' , 1.0, 'L_1' , 1.0, ...
                     'm_n2', 1.0, 'b_2' , 1.0, 'd_2' , 1.0, 'L_2' , 1.0, ...
                     'm_n3', 1.0, 'b_3' , 1.0, 'd_3' , 1.0, 'L_3' , 1.0, ...
                     'J_R' , 1.0, 'J_G' , 1.0, ...
                     'd_S' , 1.0, 'L_S' , 1.0, ...
                     'dynamic_model', @Kahraman_94};
            
            gamma = scaling_factor.process_varargin(gamma, varargin);
            dynamic_model = gamma.dynamic_model;
            gamma = rmfield(gamma, 'dynamic_model');
            gamma = scaling_factor(gamma);

            N_st = 3;
            stage = repmat(Gear_Set(), 1, N_st);
            
            for idx = 1:N_st
                stg = NREL_5MW.gear_stage(idx);
                
                gamma_idx = gamma.ends_with(num2str(idx));
                stage(idx) = stg.scale_by(gamma_idx);
            end
            
            P_r = gamma('P')* 5.0e3; % [kW], Rated power
            n_r = gamma('n')*12.1;   % [1/min.], Input speed
            
            LSS = NREL_5MW.shaft(0);
            
            inp_shaft = Shaft(LSS.d*gamma('d_S'), ...
                              LSS.L*gamma('L_S'), ...
                              LSS.bearing, ...
                              LSS.material);
            
            m_R = 110.0e3;
            J_R = 57231535.0  *gamma('J_R'); % according to [1, 3]
            m_G =     1900.0;
            J_G =      534.116*gamma('J_G');
            
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
                           'dynamic_model', dynamic_model);
%                            'dynamic_model', @Lin_Parker_99);
%                            'dynamic_model', @Dynamic_Formulation);

            obj.gamma = scaling_factor(gamma.name, gamma.value);

            [obj.S_H_val, obj.S_F_val, ...
                obj.S_shaft_val] = obj.safety_factors();
        end
    end
    
    %% Calculation:
    methods
        function [SH, SF, SShaft, calc] = safety_factor_stage(obj, idx)
            L_h    = 20*365*24;
            K_A    = 1.25;
            
            switch(idx)
                case 1
                    n_idx = [nan, nan, 0.0, obj.n_rotor];
                case 2
                    n_idx = [nan, nan, 0.0, obj.n_out(1)];
                case 3
                    n_idx = [obj.n_out(3), nan];
                otherwise
                    error();
            end
            
            calc = MATISO_6336(obj.stage(idx), 'P_rated'     , obj.P_rated, ...
                                              'n_out'       , obj.n_out(idx), ...
                                              'S_Hmin'      , obj.S_Hmin, ...
                                              'S_Fmin'      , obj.S_Fmin, ...
                                              'L_h'         , L_h, ... % [h], Required life
                                              'K_A'         , K_A, ... % [-], Application factor
                                              'lubricant_ID', 11170, ...
                                              'nu_40'       , 220.0, ... % [], lubricant viscosity
                                              'n_nominal'   , n_idx);
            [SH, SF] = calc.safety_factors('lubricant_ID', 11220, ...
                                           'nu_40'       , 220.0, ...
                                           'stage_idx'   ,   idx, ...
                                           'save_report' , false, ...
                                           'show_report' , false);
            
            K_f  = 1.0; % [-], Fatigue stress-concentration factor for bending
            K_fs = 1.0; % [-], Fatigue stress-concentration factor for torsion
            
            SShaft = obj.stage(idx).output_shaft.safety_factors(K_f, K_fs, obj.T_out(idx));
        end
    end
    
    methods(Static)
        function brg = bearing(idx)
            %BEARING returns the bearing sets of each stage of the NREL 
            % 5 MW wind turbine drivetrain according to [1], tables VI and
            % IX shown below in a combined form. Damping values were taken
            % from a SIMPACK model provided by A. Nejad. Labels were
            % obtained from manufacturers catalogs from the geometric
            % dimensions of each bearing. See remarks at the end of the
            % table.
            %
            % References: (access on 16.07.2020)
            % - http://bearingsize.info/
            % - https://simplybearings.co.uk/shop/
            % - https://www.skf.com/group/products/rolling-bearings
            %
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    Name   |           Label          | K_x (N/m) | K_y (N/m) | K_z (N/m) | K_beta (Nm/rad) | K_gamma (Nm/rad) | Type |  OD  |  ID  |  B  |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   INP-A   |     SKF C 30/1250 MB*    |    0.0    |   1.5e10  |   1.5e10  |       5e6       |        5e6       | CARB | 1750 | 1250 | 375 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   INP-B   |    SKF 231/750 CA/W33    |   4.06e8  |  1.54e10  |  1.54e10  |        0        |         0        |  SRB | 1220 |  750 | 365 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   PLC-A   |    SKF 240/710 ECA/W33   |   6.6e4   |   1.7e9   |   1.1e9   |      5.6e5      |       1.3e5      |  SRB | 1030 |  710 | 315 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   PLC-B   | SKF NF 28/1000 ECMP/HA1+ |   6.6e7   |   1.7e9   |   1.1e9   |      5.6e5      |       1.3e5      |  CRB | 1220 | 1000 | 128 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    PL-A   |     SKF NNCF 5080 CV     |   9.1e4   |   9.4e9   |   3.2e9   |      1.4e6      |       4.5e6      |  CRB |  600 |  400 | 272 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    PL-B   |     SKF NNCF 5080 CV     |   9.1e4   |   9.4e9   |   3.2e9   |      1.4e6      |       4.5e6      |  CRB |  600 |  400 | 272 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % | IMS-PLC-A |       SKF C30/710M       |   9.1e4   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      | CARB | 1030 |  710 | 236 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % | IMS-PLC-B |    SKF NU 20/600 ECMA    |   9.1e7   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  870 |  600 | 155 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |  IMS-PL-A |     SKF NNCF 4976 CV+    |   9.1e4   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  520 |  380 | 140 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |  IMS-PL-B |     SKF NNCF 4976 CV+    |   9.1e4   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  520 |  380 | 140 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   IMS-A   |     SKF NNCF 4880 CV+    |     0     |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  500 |  400 | 100 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   IMS-B   |  Timken M667948/M667911& |   7.4e7   |    5e8    |    5e8    |      1.6e6      |       1.8e6      |  TRB |  550 |  410 |  86 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   IMS-C   |  Timken M667948/M667911& |   7.8e7   |   7.4e8   |   3.3e8   |      1.1e6      |       2.5e6      |  TRB |  550 |  410 |  86 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    HS-A   |     SKF NCF 2240 ECJB    |   1.3e8   |   8.2e8   |   8.2e8   |      1.7e5      |        1e6       |  CRB |  360 |  200 |  98 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    HS-B   |        SKF 32060 X       |   6.7e7   |    8e8    |   1.3e8   |      1.7e5      |        1e6       |  TRB |  460 | 200? | 100 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    HS-C   |        SKF 32060 X       |    8e7    |    1e9    |   7.3e7   |      1.7e5      |        1e6       |  TRB |  460 | 200? | 100 |
            % +-----------+--------------------------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % *: not found on SKF website (access 16.07.2020)
            % +: multiple bearings with the same geometric parameters
            % &: in imperial units and rounded values in mm
            % ?: suspect of a typo, the correct would be 300 instead of 200
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
                    IMS_PL_A  = Bearing('name',  'IMS_PL_A', 'type'  , 'CRB' , ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'OD'  ,       520.0, 'ID'    ,  380.0, 'B'      , 140.0);
                    IMS_PL_B  = Bearing('name',  'IMS_PL_B', 'type'  , 'CRB' , ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'OD'  ,       520.0, 'ID'    ,  380.0, 'B'      , 140.0);
                    IMS_PLC_A = Bearing('name', 'IMS_PLC_A', 'type'  , 'CARB', ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'OD'  ,      1030.0, 'ID'    ,  710.0, 'B'      , 236.0);
                    IMS_PLC_B = Bearing('name', 'IMS_PLC_B', 'type'  , 'CRB' , ...
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
                    IMS_A     = Bearing('name', 'IMS_A', 'type'  , 'CRB', ...
                                        'K_x' ,     0.0, 'K_y'   , 6.0e7, 'K_z'    , 1.2e9, ...
                                                         'K_beta', 7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   360.0, 'ID'    , 200.0, 'B'      ,  98.0);
                    IMS_B     = Bearing('name', 'IMS_B', 'type'  , 'TRB', ...
                                        'K_x' ,   7.4e7, 'K_y'   , 5.0e8, 'K_z'    , 5.0e8, ...
                                                         'K_beta', 1.6e6, 'K_gamma', 1.8e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   460.0, 'ID'    , 200.0, 'B'      , 100.0);
                    IMS_C     = Bearing('name', 'IMS_C', 'type'  , 'TRB', ...
                                        'K_x' ,   7.8e7, 'K_y'   , 7.4e8, 'K_z'    , 3.3e8, ...
                                                         'K_beta', 1.1e6, 'K_gamma', 2.5e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   460.0, 'ID'    , 200.0, 'B'      , 100.0);
                    HS_A      = Bearing('name',  'HS_A', 'type'  , 'CRB', ...
                                        'K_x' ,   1.3e8, 'K_y'   , 8.2e8, 'K_z'    , 8.2e8, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   500.0, 'ID'    , 400.0, 'B'      , 100.0);
                    HS_B      = Bearing('name',  'HS_B', 'type'  , 'TRB', ...
                                        'K_x' ,   6.7e7, 'K_y'   , 8.0e8, 'K_z'    , 1.3e8, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   550.0, 'ID'    , 410.0, 'B'      ,  86.0);
                    HS_C      = Bearing('name',  'HS_C', 'type'  , 'TRB', ...
                                        'K_x' ,   8.0e7, 'K_y'   , 1.0e9, 'K_z'    , 7.3e7, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'OD'  ,   550.0, 'ID'    , 410.0, 'B'      ,  86.0);
                    
                    brg = [HS_A,  HS_B,  HS_C, ... % Pinion
                           IMS_A, IMS_B, IMS_C];   % Wheel

                otherwise
                    error('NREL_5MW:bearing', 'Option [%d] is NOT valid.', idx);
            end
        end
        
        function sft = shaft(idx)
            %SHAFT returns the shafts of each stage of the NREL 5 MW
            % wind turbine drivetrain according to [1]. Based mainly on the
            % SIMPACK simulation provided by its first author.
            %
            
            bearing = Bearing();
            
            switch(idx)
                case 0 % LSS
%                     d_0 = 1000.0;
%                     L_0 = 3000.0;
                    d = 700.0;
                    L = 2000.0;
                    bearing = NREL_5MW.bearing(0);
                    
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
                    error('NREL_5MW:shaft', 'Option [%d] is NOT valid.', idx);
            end
            
            sft = Shaft(d, L, bearing, Material());
        end
        
        function g_set = gear_stage(idx)
            %GEAR_STAGE returns the gear stages of the NREL 5 MW wind
            % turbine drivetrain according to [4], Table V reproduced below.
            % The values for the tip alteration coefficients were taken 
            % from KISSsoft.
            %
            % +----------------------------------------+-------------+--------------+-------------+
            % | Parameter                              | First stage | Second stage | Third stage |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Type                                   |  Planetary  |   Planetary  |   Parallel  |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Ratio                                  |   1:3.947   |    1:6.167   |   1:3.958   |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Number of planets                      |      3      |       3      |      -      |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Normal module (mm)                     |      45     |      21      |      14     |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Normal pressure angle (degree)         |      20     |      20      |      20     |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Helix angle (degree)                   |      0      |       0      |      10     |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Face width (mm)                        |     491     |      550     |     360     |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Centre distance (mm)                   |     863     |      584     |     861     |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Number of teeth, sun/pinion            |      19     |      18      |      24     |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Number of teeth, planet/gear           |      17     |      36      |      95     |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Number of teeth, ring gear             |      56     |      93      |      -      |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Profile shift coefficient, sun/pinion  |    0.617    |     0.389    |    0.480    |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Profile shift coefficient, planet/gear |    0.802    |     0.504    |    0.669    |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Profile shift coefficient, ring gear   |    -0.501   |     0.117    |      -      |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Pitch diameter (mm), sun/pinion        |   855.000   |    378.000   |   341.183   |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Pitch diameter (mm), planet/gear       |   765.000   |    756.000   |   1350.517  |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Pitch diameter (mm), ring gear         |   2520.000  |   1953.000   |      —      |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Tip diameter (mm), sun/pinion          |   978.839   |    432.845   |   380.751   |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Tip diameter (mm), planet/gear         |   905.440   |    815.655   |   1395.372  |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Tip diameter (mm), ring gear           |   2475.087  |   1906.074   |      -      |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Root diameter (mm), sun/pinion         |   798.061   |    341.845   |   319.627   |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Root diameter (mm), planet/gear        |   724.662   |    724.655   |   1334.248  |
            % +----------------------------------------+-------------+--------------+-------------+
            % | Root diameter (mm), ring gear          |   2677.507  |   2000.574   |      -      |
            % +----------------------------------------+-------------+--------------+-------------+
            %
            
            mat = repmat(Material(), 1, 3);
            
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
                    config = 'planetary';
                    
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
                    config = 'planetary';

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
                    config = 'parallel';
                    mat = mat(1:2);
                    
                otherwise
                    error('NREL_5MW:gear_stage', 'Option [%d] is NOT valid.', idx);
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
                             'shaft'        , NREL_5MW.shaft(idx), ...
                             'material'     , mat);
        end
        
        function property_estimation()
            u = 97.0;
            m_r = 110.0e3; % [kg], according to [2], Table 1-1
            J_g = 534.116; % [kg-m^2], according to [2], Table 5-1
            k_eq = 867637.0e3; % [N-m/rad], taken from [2], Table 5-1
            
            freq_free = 2.18; % [Hz], taken from [3], p. 45
            freq_fix  = mean([0.6205 0.6094]); % [Hz], taken from [2], Table 9-1. Mean of the results obtained using FAST and ADAMS, respectively.
            
            rho = Material().rho*1.0e9;
            
            J_r_fix   = k_eq/(2*pi*freq_fix)^2;
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
    
    methods
        function update_subvar_rigid(obj, varargin)
            %UPDATE_SUBVAR update some values of the .Rigid.SubVar file 
            % from the SIMPACK model of the NREL_5MW Drivetrain.
            %
            % Kinematic Tree:
            % ================================================================================
            %     $R_Isys
            %      ---->$B_bed_plate................................................... [0 dof]  Joint $J_bed_plate of type 35
            %           ---->$B_main_shaft............................................. [6 dof]  Joint $J_main_shaft of type 20
            %                ---->$BG_stage_01.$B_carrier.............................. [0 dof]  Joint $BG_stage_01.$J_carrier of type 0
            %                     ---->$BG_stage_01.$BG_planet_01.$B_pin............... [6 dof]  Joint $BG_stage_01.$BG_planet_01.$J_pin of type 20
            %                          ---->$BG_stage_01.$BG_planet_01.$B_gear......... [0 dof]  Joint $BG_stage_01.$BG_planet_01.$J_gear of type 0
            %
            %                     ---->$BG_stage_01.$BG_planet_02.$B_pin............... [6 dof]  Joint $BG_stage_01.$BG_planet_02.$J_pin of type 20
            %                          ---->$BG_stage_01.$BG_planet_02.$B_gear......... [0 dof]  Joint $BG_stage_01.$BG_planet_02.$J_gear of type 0
            %
            %                     ---->$BG_stage_01.$BG_planet_03.$B_pin............... [6 dof]  Joint $BG_stage_01.$BG_planet_03.$J_pin of type 20
            %                          ---->$BG_stage_01.$BG_planet_03.$B_gear......... [0 dof]  Joint $BG_stage_01.$BG_planet_03.$J_gear of type 0
            %
            % --        ---->$B_GB_frame............................................... [0 dof]  Joint $J_GB_frame of type 0
            % |              ---->$BG_stage_01.$B_ring................................. [0 dof]  Joint $BG_stage_01.$J_ring of type 0
            % |              ---->$BG_stage_01.$B_shaft................................ [6 dof]  Joint $BG_stage_01.$J_shaft of type 20
            % |                   ---->$BG_stage_01.$B_sun............................. [0 dof]  Joint $BG_stage_01.$J_sun of type 0
            % |                   ---->$BG_stage_02.$B_carrier......................... [0 dof]  Joint $BG_stage_02.$J_carrier of type 0
            % |                        ---->$BG_stage_02.$BG_planet_01.$B_pin.......... [6 dof]  Joint $BG_stage_02.$BG_planet_01.$J_pin of type 20
            % |                             ---->$BG_stage_02.$BG_planet_01.$B_gear.... [0 dof]  Joint $BG_stage_02.$BG_planet_01.$J_gear of type 0
            % |
            % |                        ---->$BG_stage_02.$BG_planet_02.$B_pin.......... [6 dof]  Joint $BG_stage_02.$BG_planet_02.$J_pin of type 20
            % |                             ---->$BG_stage_02.$BG_planet_02.$B_gear.... [0 dof]  Joint $BG_stage_02.$BG_planet_02.$J_gear of type 0
            % |
            % |                        ---->$BG_stage_02.$BG_planet_03.$B_pin.......... [6 dof]  Joint $BG_stage_02.$BG_planet_03.$J_pin of type 20
            % |                             ---->$BG_stage_02.$BG_planet_03.$B_gear.... [0 dof]  Joint $BG_stage_02.$BG_planet_03.$J_gear of type 0
            % |
            % |              ---->$BG_stage_02.$B_ring................................. [0 dof]  Joint $BG_stage_02.$J_ring of type 0
            % |              ---->$BG_stage_02.$B_shaft................................ [6 dof]  Joint $BG_stage_02.$J_shaft of type 20
            % |                   ---->$BG_stage_02.$B_sun............................. [0 dof]  Joint $BG_stage_02.$J_sun of type 0
            % |                   ---->$BG_stage_03.$BG_wheel.$B_pin................... [0 dof]  Joint $BG_stage_03.$BG_wheel.$J_pin of type 0
            % |                        ---->$BG_stage_03.$BG_wheel.$B_gear............. [0 dof]  Joint $BG_stage_03.$BG_wheel.$J_gear of type 0
            % |
            % --             ---->$BG_stage_03.$B_shaft................................ [6 dof]  Joint $BG_stage_03.$J_shaft of type 20
            % .
            % $L_output_shaft.......................................................... [1 cft]  Constraint of type 25
            %                     ---->$BG_stage_03.$BG_pinion.$B_pin.................. [0 dof]  Joint $BG_stage_03.$BG_pinion.$J_pin of type 0
            %                          ---->$BG_stage_03.$BG_pinion.$B_gear............ [0 dof]  Joint $BG_stage_03.$BG_pinion.$J_gear of type 0
            %
            
            if(isempty(varargin{:}))
                mesh_flag = false;
            else
                mesh_flag = varargin{1};
            end
            
            gamma_x = max([obj.gamma('m_n1') obj.gamma('m_n2') obj.gamma('m_n3')]);
            gamma_load = gamma_x^2;
            
            k_sp = zeros(2, 1);
            k_rp = k_sp;
            for idx = 1:2
                k_sp(idx) = obj.stage(idx).sub_set('sun_planet').k_mesh;
                k_rp(idx) = obj.stage(idx).sub_set('planet_ring').k_mesh;
            end
            
            k_3 = obj.stage(idx + 1).k_mesh;
            
            new_ID = fopen('@NREL_5MW\NREL_5MW.rigid.subvar', 'w');
            
            fprintf(new_ID, "!file.version=3.4! Removing this line will make the file unreadable\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "!**********************************************************************\n");
            fprintf(new_ID, "! SubVars\n");
            fprintf(new_ID, "!**********************************************************************\n");
            fprintf(new_ID, "subvar($_rated_speed, str = '%e rpm')                                                        ! $_rated_speed\n", obj.n_rotor);
            fprintf(new_ID, "subvar($_damping_coeff, str = '10e-6 m')                                                       ! $_damping_coeff\n");
            fprintf(new_ID, "subvar($_addendum_coeff, str = '1.0')                                                          ! $_addendum_coeff\n");
            fprintf(new_ID, "subvar($_dedendum_coeff, str = '1.25')                                                         ! $_dedendum_coeff\n");
            fprintf(new_ID, "subvar($_regularization_vel, str = '1.0e-3')                                                   ! $_regularization_vel\n");
            fprintf(new_ID, "subvar($_tooth_damping, str = '5.0e8 N/(m/s)')                                                 ! $_tooth_damping\n");
            fprintf(new_ID, "subvar($_tooth_stiff_ratio, str = '0.8')                                                       ! $_tooth_stiff_ratio\n");
            fprintf(new_ID, "subvar($_height, str = '1.75 * %e m')                                                               ! $_height\n", gamma_x);
            fprintf(new_ID, "subvar($_gamma_load, str = '%e')                                                                ! $_gamma_load\n", gamma_load);
            fprintf(new_ID, "subvar($_gamma_x, str = '%e')                                                                   ! $_gamma_x\n", gamma_x);
            fprintf(new_ID, "subvargroup.begin (                 $SVG_material                 )                           ! $SVG_material\n");
            fprintf(new_ID, "   subvar($_E, str = '206.0e9 Pa', desc (   1) = 'Young\\'s modulus')                          ! $SVG_material.$_E\n");
            fprintf(new_ID, "   subvar($_nu, str = '0.3', desc (   1) = 'Poisson\\'s ratio')                                ! $SVG_material.$_nu\n");
            fprintf(new_ID, "   subvar($_rho, str = '7.83e3 kg/m^3', desc (   1) = 'Density')                               ! $SVG_material.$_rho\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_material                 )                           ! $SVG_material\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_bed_plate                )                           ! $SVG_bed_plate\n");
            fprintf(new_ID, "   subvar($_length, str = 'anint($SVG_bed_plate.$_x_gen) m')                                   ! $SVG_bed_plate.$_length\n");
            fprintf(new_ID, "   subvar($_width, str = '$SVG_bed_plate.$_length')                                            ! $SVG_bed_plate.$_width\n");
            fprintf(new_ID, "   subvar($_thickness, str = '0.5 * %e m')                                                          ! $SVG_bed_plate.$_thickness\n", gamma_x);
            fprintf(new_ID, "   subvar($_x_gen, str = '$SVG_stage_03.$SVG_wheel.$SVG_C.$_x_GB + $SVG_stage_03.$SVG_shaft.$_length') ! $SVG_bed_plate.$_x_gen\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_bed_plate                )                           ! $SVG_bed_plate\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_main_shaft               )                           ! $SVG_main_shaft\n");
            fprintf(new_ID, "   subvar($_diameter, str = '%e mm')                                                        ! $SVG_main_shaft.$_diameter\n", obj.main_shaft.d);
            fprintf(new_ID, "   subvar($_length, str = '%e mm')                                                         ! $SVG_main_shaft.$_length\n", obj.main_shaft.L);
            fprintf(new_ID, "subvargroup.end (                   $SVG_main_shaft               )                           ! $SVG_main_shaft\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_GB_frame                 )                           ! $SVG_GB_frame\n");
            fprintf(new_ID, "   subvar($_length, str = '3.5 * %e m')                                                             ! $SVG_GB_frame.$_length\n", gamma_x);
            fprintf(new_ID, "   subvar($_width, str = '2.0*$_height')                                                       ! $SVG_GB_frame.$_width\n");
            fprintf(new_ID, "   subvar($_x, str = '($SVG_main_shaft.$_length*7/4 - $SVG_bed_plate.$_length + $SVG_GB_frame.$_length)/2.0') ! $SVG_GB_frame.$_x\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_GB_frame                 )                           ! $SVG_GB_frame\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_INP_A                    )                           ! $SVG_INP_A\n");
            fprintf(new_ID, "   subvar($_x, str = '$SVG_main_shaft.$_length/6.0')                                           ! $SVG_INP_A.$_x\n");
            fprintf(new_ID, "   subvar($_K_x, str = '0.0 N/m')                                                              ! $SVG_INP_A.$_K_x\n");
            fprintf(new_ID, "   subvar($_K_y, str = '%e N/m')                                                          ! $SVG_INP_A.$_K_y\n", obj.main_shaft.bearing(1).K_y);
            fprintf(new_ID, "   subvar($_K_z, str = '$SVG_INP_A.$_K_y')                                                     ! $SVG_INP_A.$_K_z\n");
            fprintf(new_ID, "   subvar($_K_alpha, str = '0.0 Nm/rad')                                                       ! $SVG_INP_A.$_K_alpha\n");
            fprintf(new_ID, "   subvar($_K_beta, str = '%e Nm/rad')                                                      ! $SVG_INP_A.$_K_beta\n", obj.main_shaft.bearing(1).K_beta);
            fprintf(new_ID, "   subvar($_K_gamma, str = '$SVG_INP_A.$_K_beta')                                              ! $SVG_INP_A.$_K_gamma\n");
            fprintf(new_ID, "   subvar($_D_x, str = '0.0 N/(m/s)')                                                          ! $SVG_INP_A.$_D_x\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_INP_A                    )                           ! $SVG_INP_A\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_INP_B                    )                           ! $SVG_INP_B\n");
            fprintf(new_ID, "   subvar($_x, str = '$SVG_main_shaft.$_length*3.0/4.0')                                       ! $SVG_INP_B.$_x\n");
            fprintf(new_ID, "   subvar($_K_x, str = '%e N/m')                                                           ! $SVG_INP_B.$_K_x\n", obj.main_shaft.bearing(2).K_x);
            fprintf(new_ID, "   subvar($_K_y, str = '%e N/m')                                                          ! $SVG_INP_B.$_K_y\n", obj.main_shaft.bearing(2).K_y);
            fprintf(new_ID, "   subvar($_K_z, str = '$SVG_INP_B.$_K_y')                                                     ! $SVG_INP_B.$_K_z\n");
            fprintf(new_ID, "   subvar($_K_alpha, str = '0.0 Nm/rad')                                                       ! $SVG_INP_B.$_K_alpha\n");
            fprintf(new_ID, "   subvar($_K_beta, str = '0.0 Nm/rad')                                                        ! $SVG_INP_B.$_K_beta\n");
            fprintf(new_ID, "   subvar($_K_gamma, str = '$SVG_INP_B.$_K_beta')                                              ! $SVG_INP_B.$_K_gamma\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_INP_B                    )                           ! $SVG_INP_B\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_stage_01                 )                           ! $SVG_stage_01\n");
            fprintf(new_ID, "   subvar($_ring_planet_stiffness, str = '%e N/m')                                          ! $SVG_stage_01.$_ring_planet_stiffness\n", k_rp(1));
            fprintf(new_ID, "   subvar($_sun_planet_stiffness, str = '%e N/m')                                           ! $SVG_stage_01.$_sun_planet_stiffness\n", k_sp(1));
            fprintf(new_ID, "   subvar($_normal_module, str = '%e mm')                                                    ! $SVG_stage_01.$_normal_module\n", obj.stage(1).m_n);
            fprintf(new_ID, "   subvar($_face_width, str = '%e mm')                                                      ! $SVG_stage_01.$_face_width\n", obj.stage(1).b);
            fprintf(new_ID, "   subvar($_center_distance, str = '%e mm')                                                 ! $SVG_stage_01.$_center_distance\n", obj.stage(1).a_w);
            fprintf(new_ID, "   subvar($_pressure_angle, str = '20.0 deg')                                                  ! $SVG_stage_01.$_pressure_angle\n");
            fprintf(new_ID, "   subvar($_helix_angle, str = '0.0 deg')                                                      ! $SVG_stage_01.$_helix_angle\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_sun                      )                           ! $SVG_stage_01.$SVG_sun\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '19')                                                          ! $SVG_stage_01.$SVG_sun.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '0.617')                                                ! $SVG_stage_01.$SVG_sun.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_01.$SVG_sun.$_diam_bore\n", obj.stage(1).d_b(1) * ...
                obj.stage(1).bore_ratio(1));
            fprintf(new_ID, "      subvar($_angle, str = '11.5*(180deg)/$SVG_stage_01.$SVG_sun.$_num_teeth')                ! $SVG_stage_01.$SVG_sun.$_angle\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_sun                      )                           ! $SVG_stage_01.$SVG_sun\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_planet                   )                           ! $SVG_stage_01.$SVG_planet\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '17')                                                          ! $SVG_stage_01.$SVG_planet.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '0.802')                                                ! $SVG_stage_01.$SVG_planet.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_01.$SVG_planet.$_diam_bore\n", obj.stage(1).d_b(2) * ...
                obj.stage(1).bore_ratio(2));
            fprintf(new_ID, "      subvar($_pin_width, str = '1.4*$SVG_stage_01.$_face_width')                              ! $SVG_stage_01.$SVG_planet.$_pin_width\n");
            fprintf(new_ID, "      subvar($_angle, str = '5.25*(360deg)/$SVG_stage_01.$SVG_planet.$_num_teeth')             ! $SVG_stage_01.$SVG_planet.$_angle\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_A                        )                           ! $SVG_stage_01.$SVG_planet.$SVG_A\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_01.$SVG_carrier.$_width/2.0')                           ! $SVG_stage_01.$SVG_planet.$SVG_A.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_01.$SVG_planet.$SVG_A.$_K_x\n", obj.stage(1).bearing(1).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_01.$SVG_planet.$SVG_A.$_K_y\n", obj.stage(1).bearing(1).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_01.$SVG_planet.$SVG_A.$_K_z\n", obj.stage(1).bearing(1).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_01.$SVG_planet.$SVG_A.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_01.$SVG_planet.$SVG_A.$_K_beta\n", obj.stage(1).bearing(1).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '%e Nm/rad')                                               ! $SVG_stage_01.$SVG_planet.$SVG_A.$_K_gamma\n", obj.stage(1).bearing(1).K_gamma);
            fprintf(new_ID, "      subvargroup.end (             $SVG_A                        )                           ! $SVG_stage_01.$SVG_planet.$SVG_A\n");
            fprintf(new_ID, "      subvargroup.desc (   1 ,      $SVG_A                        ) = 'PL_A, PL_B, and PL_C'  ! Description\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_planet                   )                           ! $SVG_stage_01.$SVG_planet\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_ring                     )                           ! $SVG_stage_01.$SVG_ring\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '56')                                                          ! $SVG_stage_01.$SVG_ring.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '-0.501')                                               ! $SVG_stage_01.$SVG_ring.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_01.$SVG_ring.$_diam_bore\n", obj.stage(1).d_b(3) * ...
                obj.stage(1).bore_ratio(3));
            fprintf(new_ID, "      subvar($_angle, str = '27*(180 deg)/$SVG_stage_01.$SVG_ring.$_num_teeth')                ! $SVG_stage_01.$SVG_ring.$_angle\n");
            fprintf(new_ID, "      subvar($_x, str = '$SVG_main_shaft.$_length + $SVG_stage_01.$SVG_carrier.$_width/2.0')   ! $SVG_stage_01.$SVG_ring.$_x\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_ring                     )                           ! $SVG_stage_01.$SVG_ring\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_carrier                  )                           ! $SVG_stage_01.$SVG_carrier\n");
            fprintf(new_ID, "      subvar($_diameter_in, str = '1.4*$SVG_stage_01.$_center_distance	')                      ! $SVG_stage_01.$SVG_carrier.$_diameter_in\n");
            fprintf(new_ID, "      subvar($_diameter_out, str = '2.6*$SVG_stage_01.$_center_distance	')                     ! $SVG_stage_01.$SVG_carrier.$_diameter_out\n");
            fprintf(new_ID, "      subvar($_width, str = '1.2*$SVG_stage_01.$_face_width')                                  ! $SVG_stage_01.$SVG_carrier.$_width\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_A                        )                           ! $SVG_stage_01.$SVG_carrier.$SVG_A\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_01.$SVG_carrier.$_width/2.0')                           ! $SVG_stage_01.$SVG_carrier.$SVG_A.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_01.$SVG_carrier.$SVG_A.$_K_x\n", obj.stage(1).bearing(3).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_01.$SVG_carrier.$SVG_A.$_K_y\n", obj.stage(1).bearing(3).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_01.$SVG_carrier.$SVG_A.$_K_z\n", obj.stage(1).bearing(3).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_01.$SVG_carrier.$SVG_A.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_01.$SVG_carrier.$SVG_A.$_K_beta\n", obj.stage(1).bearing(3).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '%e Nm/rad')                                               ! $SVG_stage_01.$SVG_carrier.$SVG_A.$_K_gamma\n", obj.stage(1).bearing(3).K_gamma);
            fprintf(new_ID, "      subvargroup.end (             $SVG_A                        )                           ! $SVG_stage_01.$SVG_carrier.$SVG_A\n");
            fprintf(new_ID, "      subvargroup.desc (   1 ,      $SVG_A                        ) = 'PLC_A'                 ! Description\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_B                        )                           ! $SVG_stage_01.$SVG_carrier.$SVG_B\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_main_shaft.$_length+$SVG_stage_01.$SVG_carrier.$_width')      ! $SVG_stage_01.$SVG_carrier.$SVG_B.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_01.$SVG_carrier.$SVG_B.$_K_x\n", obj.stage(1).bearing(4).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '$SVG_stage_01.$SVG_carrier.$SVG_A.$_K_y')                        ! $SVG_stage_01.$SVG_carrier.$SVG_B.$_K_y\n");
            fprintf(new_ID, "         subvar($_K_z, str = '$SVG_stage_01.$SVG_carrier.$SVG_A.$_K_z')                        ! $SVG_stage_01.$SVG_carrier.$SVG_B.$_K_z\n");
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_01.$SVG_carrier.$SVG_B.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '$SVG_stage_01.$SVG_carrier.$SVG_A.$_K_beta')                  ! $SVG_stage_01.$SVG_carrier.$SVG_B.$_K_beta\n");
            fprintf(new_ID, "         subvar($_K_gamma, str = '$SVG_stage_01.$SVG_carrier.$SVG_A.$_K_gamma')                ! $SVG_stage_01.$SVG_carrier.$SVG_B.$_K_gamma\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_B                        )                           ! $SVG_stage_01.$SVG_carrier.$SVG_B\n");
            fprintf(new_ID, "      subvargroup.desc (   1 ,      $SVG_B                        ) = 'PLC_B'                 ! Description\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_carrier                  )                           ! $SVG_stage_01.$SVG_carrier\n");
            fprintf(new_ID, "   subvargroup.desc (     1 ,       $SVG_carrier                  ) = 'PLC_A and PLC_B'       ! Description\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_shaft                    )                           ! $SVG_stage_01.$SVG_shaft\n");
            fprintf(new_ID, "      subvar($_diameter, str = '%e mm')                                                     ! $SVG_stage_01.$SVG_shaft.$_diameter\n", obj.stage(1).output_shaft.d);
            fprintf(new_ID, "      subvar($_length, str = '%e mm')                                                       ! $SVG_stage_01.$SVG_shaft.$_length\n", obj.stage(1).output_shaft.L);
            fprintf(new_ID, "   subvargroup.end (                $SVG_shaft                    )                           ! $SVG_stage_01.$SVG_shaft\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_stage_01                 )                           ! $SVG_stage_01\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_stage_02                 )                           ! $SVG_stage_02\n");
            fprintf(new_ID, "   subvar($_ring_planet_stiffness, str = '%e N/m')                                          ! $SVG_stage_01.$_ring_planet_stiffness\n", k_rp(2));
            fprintf(new_ID, "   subvar($_sun_planet_stiffness, str = '%e N/m')                                           ! $SVG_stage_01.$_sun_planet_stiffness\n", k_sp(2));
            fprintf(new_ID, "   subvar($_normal_module, str = '%e mm')                                                    ! $SVG_stage_02.$_normal_module\n", obj.stage(2).m_n);
            fprintf(new_ID, "   subvar($_face_width, str = '%e mm')                                                      ! $SVG_stage_02.$_face_width\n", obj.stage(2).b);
            fprintf(new_ID, "   subvar($_center_distance, str = '%e mm')                                                 ! $SVG_stage_02.$_center_distance\n", obj.stage(2).a_w);
            fprintf(new_ID, "   subvar($_pressure_angle, str = '20.0 deg')                                                  ! $SVG_stage_02.$_pressure_angle\n");
            fprintf(new_ID, "   subvar($_helix_angle, str = '0.0 deg')                                                      ! $SVG_stage_02.$_helix_angle\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_sun                      )                           ! $SVG_stage_02.$SVG_sun\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '18')                                                          ! $SVG_stage_02.$SVG_sun.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '0.389')                                                ! $SVG_stage_02.$SVG_sun.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_02.$SVG_sun.$_diam_bore\n", obj.stage(2).d_b(1) * ...
                obj.stage(2).bore_ratio(1));
            fprintf(new_ID, "      subvar($_angle, str = '12*(180deg)/$SVG_stage_02.$SVG_sun.$_num_teeth')                  ! $SVG_stage_02.$SVG_sun.$_angle\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_sun                      )                           ! $SVG_stage_02.$SVG_sun\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_planet                   )                           ! $SVG_stage_02.$SVG_planet\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '36')                                                          ! $SVG_stage_02.$SVG_planet.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '0.504')                                                ! $SVG_stage_02.$SVG_planet.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_02.$SVG_planet.$_diam_bore\n", obj.stage(2).d_b(2) * ...
                obj.stage(2).bore_ratio(2));
            fprintf(new_ID, "      subvar($_pin_width, str = '1.4*$SVG_stage_02.$_face_width')                              ! $SVG_stage_02.$SVG_planet.$_pin_width\n");
            fprintf(new_ID, "      subvar($_angle, str = '20*(180 deg)/$SVG_stage_02.$SVG_planet.$_num_teeth')              ! $SVG_stage_02.$SVG_planet.$_angle\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_A                        )                           ! $SVG_stage_02.$SVG_planet.$SVG_A\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_02.$SVG_carrier.$_width/2.0')                           ! $SVG_stage_02.$SVG_planet.$SVG_A.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_02.$SVG_planet.$SVG_A.$_K_x\n", obj.stage(2).bearing(1).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_02.$SVG_planet.$SVG_A.$_K_y\n", obj.stage(2).bearing(1).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_02.$SVG_planet.$SVG_A.$_K_z\n", obj.stage(2).bearing(1).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_02.$SVG_planet.$SVG_A.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_02.$SVG_planet.$SVG_A.$_K_beta\n", obj.stage(2).bearing(1).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '%e Nm/rad')                                               ! $SVG_stage_02.$SVG_planet.$SVG_A.$_K_gamma\n", obj.stage(2).bearing(1).K_gamma);
            fprintf(new_ID, "      subvargroup.end (             $SVG_A                        )                           ! $SVG_stage_02.$SVG_planet.$SVG_A\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_planet                   )                           ! $SVG_stage_02.$SVG_planet\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_ring                     )                           ! $SVG_stage_02.$SVG_ring\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '93')                                                          ! $SVG_stage_02.$SVG_ring.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '0.117')                                                ! $SVG_stage_02.$SVG_ring.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_02.$SVG_ring.$_diam_bore\n", obj.stage(2).d_b(3) * ...
                obj.stage(2).bore_ratio(3));
            fprintf(new_ID, "      subvar($_angle, str = '45.5*(180 deg)/$SVG_stage_02.$SVG_ring.$_num_teeth')              ! $SVG_stage_02.$SVG_ring.$_angle\n");
            fprintf(new_ID, "      subvar($_x, str = '$SVG_main_shaft.$_length + $SVG_stage_01.$SVG_carrier.$_width + $SVG_stage_01.$SVG_shaft.$_length + $SVG_stage_02.$SVG_carrier.$_width/2.0') ! $SVG_stage_02.$SVG_ring.$_x\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_ring                     )                           ! $SVG_stage_02.$SVG_ring\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_carrier                  )                           ! $SVG_stage_02.$SVG_carrier\n");
            fprintf(new_ID, "      subvar($_diameter_in, str = '1.4*$SVG_stage_02.$_center_distance	')                      ! $SVG_stage_02.$SVG_carrier.$_diameter_in\n");
            fprintf(new_ID, "      subvar($_diameter_out, str = '2.6*$SVG_stage_02.$_center_distance	')                     ! $SVG_stage_02.$SVG_carrier.$_diameter_out\n");
            fprintf(new_ID, "      subvar($_width, str = '1.2*$SVG_stage_02.$_face_width')                                  ! $SVG_stage_02.$SVG_carrier.$_width\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_A                        )                           ! $SVG_stage_02.$SVG_carrier.$SVG_A\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_main_shaft.$_length + $SVG_stage_01.$SVG_carrier.$_width + $SVG_stage_01.$SVG_shaft.$_length') ! $SVG_stage_02.$SVG_carrier.$SVG_A.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_02.$SVG_carrier.$SVG_A.$_K_x\n", obj.stage(2).bearing(3).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_02.$SVG_carrier.$SVG_A.$_K_y\n", obj.stage(2).bearing(3).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_02.$SVG_carrier.$SVG_A.$_K_z\n", obj.stage(2).bearing(3).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_02.$SVG_carrier.$SVG_A.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_02.$SVG_carrier.$SVG_A.$_K_beta\n", obj.stage(2).bearing(3).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '$SVG_stage_02.$SVG_carrier.$SVG_A.$_K_beta')                 ! $SVG_stage_02.$SVG_carrier.$SVG_A.$_K_gamma\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_A                        )                           ! $SVG_stage_02.$SVG_carrier.$SVG_A\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_B                        )                           ! $SVG_stage_02.$SVG_carrier.$SVG_B\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_02.$SVG_carrier.$SVG_A.$_x + $SVG_stage_02.$SVG_carrier.$_width') ! $SVG_stage_02.$SVG_carrier.$SVG_B.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_02.$SVG_carrier.$SVG_B.$_K_x\n", obj.stage(2).bearing(4).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '$SVG_stage_02.$SVG_carrier.$SVG_A.$_K_y')                        ! $SVG_stage_02.$SVG_carrier.$SVG_B.$_K_y\n");
            fprintf(new_ID, "         subvar($_K_z, str = '$SVG_stage_02.$SVG_carrier.$SVG_A.$_K_z')                        ! $SVG_stage_02.$SVG_carrier.$SVG_B.$_K_z\n");
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_02.$SVG_carrier.$SVG_B.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '$SVG_stage_02.$SVG_carrier.$SVG_A.$_K_beta')                  ! $SVG_stage_02.$SVG_carrier.$SVG_B.$_K_beta\n");
            fprintf(new_ID, "         subvar($_K_gamma, str = '$SVG_stage_02.$SVG_carrier.$SVG_A.$_K_beta')                 ! $SVG_stage_02.$SVG_carrier.$SVG_B.$_K_gamma\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_B                        )                           ! $SVG_stage_02.$SVG_carrier.$SVG_B\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_carrier                  )                           ! $SVG_stage_02.$SVG_carrier\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_shaft                    )                           ! $SVG_stage_02.$SVG_shaft\n");
            fprintf(new_ID, "      subvar($_diameter, str = '%e mm')                                                     ! $SVG_stage_02.$SVG_shaft.$_diameter\n", obj.stage(2).output_shaft.d);
            fprintf(new_ID, "      subvar($_length, str = '%e mm')                                                       ! $SVG_stage_02.$SVG_shaft.$_length\n", obj.stage(2).output_shaft.L);
            fprintf(new_ID, "   subvargroup.end (                $SVG_shaft                    )                           ! $SVG_stage_02.$SVG_shaft\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_stage_02                 )                           ! $SVG_stage_02\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_stage_03                 )                           ! $SVG_stage_03\n");
            fprintf(new_ID, "   subvar($_mesh_stiffness, str = '%e N/m')                                                 ! $SVG_stage_03.$_mesh_stiffness\n", k_3);
            fprintf(new_ID, "   subvar($_normal_module, str = '%e mm')                                                    ! $SVG_stage_03.$_normal_module\n", obj.stage(3).m_n);
            fprintf(new_ID, "   subvar($_face_width, str = '%e mm')                                                      ! $SVG_stage_03.$_face_width\n", obj.stage(3).b);
            fprintf(new_ID, "   subvar($_center_distance, str = '%e mm')                                                 ! $SVG_stage_03.$_center_distance\n", obj.stage(3).a_w);
            fprintf(new_ID, "   subvar($_pressure_angle, str = '20.0 deg')                                                  ! $SVG_stage_03.$_pressure_angle\n");
            fprintf(new_ID, "   subvar($_helix_angle, str = '10.0 deg')                                                     ! $SVG_stage_03.$_helix_angle\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_pinion                   )                           ! $SVG_stage_03.$SVG_pinion\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '24')                                                          ! $SVG_stage_03.$SVG_pinion.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '0.480')                                                ! $SVG_stage_03.$SVG_pinion.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_03.$SVG_pinion.$_diam_bore\n", obj.stage(3).d_b(1) * ...
                obj.stage(3).bore_ratio(1));
            fprintf(new_ID, "      subvar($_pin_width, str = '1.4*$SVG_stage_03.$_face_width')                              ! $SVG_stage_03.$SVG_pinion.$_pin_width\n");
            fprintf(new_ID, "      subvar($_angle, str = '12.5*(180.0 deg)/$SVG_stage_03.$SVG_pinion.$_num_teeth')          ! $SVG_stage_03.$SVG_pinion.$_angle\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_A                        )                           ! $SVG_stage_03.$SVG_pinion.$SVG_A\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_01.$SVG_carrier.$_width/2.0')                           ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_K_x\n", obj.stage(3).bearing(1).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_K_y\n", obj.stage(3).bearing(1).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '$SVG_stage_03.$SVG_pinion.$SVG_A.$_K_y')                         ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_K_z\n");
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_K_beta\n", obj.stage(3).bearing(1).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '%e Nm/rad')                                               ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_K_gamma\n", obj.stage(3).bearing(1).K_gamma);
            fprintf(new_ID, "         subvar($_y, str = '-$SVG_stage_03.$_center_distance + $SVG_bed_plate.$_width/2.0')    ! $SVG_stage_03.$SVG_pinion.$SVG_A.$_y\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_A                        )                           ! $SVG_stage_03.$SVG_pinion.$SVG_A\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_B                        )                           ! $SVG_stage_03.$SVG_pinion.$SVG_B\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_01.$SVG_carrier.$_width/2.0')                           ! $SVG_stage_03.$SVG_pinion.$SVG_B.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_B.$_K_x\n", obj.stage(3).bearing(2).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_B.$_K_y\n", obj.stage(3).bearing(2).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_B.$_K_z\n", obj.stage(3).bearing(2).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_03.$SVG_pinion.$SVG_B.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '$SVG_stage_03.$SVG_pinion.$SVG_A.$_K_beta')                   ! $SVG_stage_03.$SVG_pinion.$SVG_B.$_K_beta\n");
            fprintf(new_ID, "         subvar($_K_gamma, str = '$SVG_stage_03.$SVG_pinion.$SVG_A.$_K_gamma')                 ! $SVG_stage_03.$SVG_pinion.$SVG_B.$_K_gamma\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_B                        )                           ! $SVG_stage_03.$SVG_pinion.$SVG_B\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_C                        )                           ! $SVG_stage_03.$SVG_pinion.$SVG_C\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_01.$SVG_carrier.$_width/2.0')                           ! $SVG_stage_03.$SVG_pinion.$SVG_C.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_C.$_K_x\n", obj.stage(3).bearing(3).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_C.$_K_y\n", obj.stage(3).bearing(3).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_pinion.$SVG_C.$_K_z\n", obj.stage(3).bearing(3).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_03.$SVG_pinion.$SVG_C.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '$SVG_stage_03.$SVG_pinion.$SVG_A.$_K_beta')                   ! $SVG_stage_03.$SVG_pinion.$SVG_C.$_K_beta\n");
            fprintf(new_ID, "         subvar($_K_gamma, str = '$SVG_stage_03.$SVG_pinion.$SVG_A.$_K_gamma')                 ! $SVG_stage_03.$SVG_pinion.$SVG_C.$_K_gamma\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_C                        )                           ! $SVG_stage_03.$SVG_pinion.$SVG_C\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_pinion                   )                           ! $SVG_stage_03.$SVG_pinion\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_wheel                    )                           ! $SVG_stage_03.$SVG_wheel\n");
            fprintf(new_ID, "      subvar($_num_teeth, str = '95')                                                          ! $SVG_stage_03.$SVG_wheel.$_num_teeth\n");
            fprintf(new_ID, "      subvar($_prof_shift_coeff, str = '0.669')                                                ! $SVG_stage_03.$SVG_wheel.$_prof_shift_coeff\n");
            fprintf(new_ID, "      subvar($_diam_bore, str = '%e mm')                                                    ! $SVG_stage_03.$SVG_wheel.$_diam_bore\n", obj.stage(3).d_b(2) * ...
                obj.stage(3).bore_ratio(2));
            fprintf(new_ID, "      subvar($_pin_width, str = '1.4*$SVG_stage_03.$_face_width')                              ! $SVG_stage_03.$SVG_wheel.$_pin_width\n");
            fprintf(new_ID, "      subvar($_angle, str = '51*(180 deg)/$SVG_stage_03.$SVG_wheel.$_num_teeth')               ! $SVG_stage_03.$SVG_wheel.$_angle\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_A                        )                           ! $SVG_stage_03.$SVG_wheel.$SVG_A\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_03.$SVG_wheel.$_pin_width/2.0')                         ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '0.0 N/m')                                                        ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_K_x\n");
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_K_y\n", obj.stage(3).bearing(4).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_K_z\n", obj.stage(3).bearing(4).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_K_beta\n", obj.stage(3).bearing(4).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '$SVG_stage_03.$SVG_wheel.$SVG_A.$_K_beta')                   ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_K_gamma\n");
            fprintf(new_ID, "         subvar($_x_GB, str = '$SVG_stage_02.$SVG_carrier.$SVG_B.$_x + $SVG_stage_02.$SVG_shaft.$_length') ! $SVG_stage_03.$SVG_wheel.$SVG_A.$_x_GB\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_A                        )                           ! $SVG_stage_03.$SVG_wheel.$SVG_A\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_B                        )                           ! $SVG_stage_03.$SVG_wheel.$SVG_B\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_03.$SVG_wheel.$SVG_A.$_x*6.0/7.0')                      ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_K_x\n", obj.stage(3).bearing(5).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_K_y\n", obj.stage(3).bearing(5).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '$SVG_stage_03.$SVG_wheel.$SVG_B.$_K_y')                          ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_K_z\n");
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_K_beta\n", obj.stage(3).bearing(5).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '%e Nm/rad')                                               ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_K_gamma\n", obj.stage(3).bearing(5).K_gamma);
            fprintf(new_ID, "         subvar($_x_GB, str = '$SVG_stage_03.$SVG_wheel.$SVG_A.$_x_GB + $SVG_stage_03.$SVG_wheel.$SVG_A.$_x + $SVG_stage_03.$SVG_wheel.$SVG_B.$_x') ! $SVG_stage_03.$SVG_wheel.$SVG_B.$_x_GB\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_B                        )                           ! $SVG_stage_03.$SVG_wheel.$SVG_B\n");
            fprintf(new_ID, "      subvargroup.begin (           $SVG_C                        )                           ! $SVG_stage_03.$SVG_wheel.$SVG_C\n");
            fprintf(new_ID, "         subvar($_x, str = '$SVG_stage_03.$SVG_wheel.$SVG_A.$_x*6.0/7.0')                      ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_x\n");
            fprintf(new_ID, "         subvar($_K_x, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_K_x\n", obj.stage(3).bearing(6).K_x);
            fprintf(new_ID, "         subvar($_K_y, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_K_y\n", obj.stage(3).bearing(6).K_y);
            fprintf(new_ID, "         subvar($_K_z, str = '%e N/m')                                                      ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_K_z\n", obj.stage(3).bearing(6).K_z);
            fprintf(new_ID, "         subvar($_K_alpha, str = '0.0 Nm/rad')                                                 ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_K_alpha\n");
            fprintf(new_ID, "         subvar($_K_beta, str = '%e Nm/rad')                                                ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_K_beta\n", obj.stage(3).bearing(6).K_beta);
            fprintf(new_ID, "         subvar($_K_gamma, str = '%e Nm/rad')                                               ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_K_gamma\n", obj.stage(3).bearing(6).K_gamma);
            fprintf(new_ID, "         subvar($_x_GB, str = '$SVG_stage_03.$SVG_wheel.$SVG_A.$_x_GB + $SVG_stage_03.$SVG_wheel.$_pin_width') ! $SVG_stage_03.$SVG_wheel.$SVG_C.$_x_GB\n");
            fprintf(new_ID, "      subvargroup.end (             $SVG_C                        )                           ! $SVG_stage_03.$SVG_wheel.$SVG_C\n");
            fprintf(new_ID, "   subvargroup.end (                $SVG_wheel                    )                           ! $SVG_stage_03.$SVG_wheel\n");
            fprintf(new_ID, "   subvargroup.begin (              $SVG_shaft                    )                           ! $SVG_stage_03.$SVG_shaft\n");
            fprintf(new_ID, "      subvar($_diameter, str = '%e mm')                                                     ! $SVG_stage_03.$SVG_shaft.$_diameter\n", obj.stage(3).output_shaft.d);
            fprintf(new_ID, "      subvar($_length, str = '%e mm')                                                      ! $SVG_stage_03.$SVG_shaft.$_length\n", obj.stage(3).output_shaft.L);
            fprintf(new_ID, "   subvargroup.end (                $SVG_shaft                    )                           ! $SVG_stage_03.$SVG_shaft\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_stage_03                 )                           ! $SVG_stage_03\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_bearing_damping          )                           ! $SVG_bearing_damping\n");
            fprintf(new_ID, "   subvar($_D_x, str = '%e N/(m/s)')                                                        ! $SVG_bearing_damping.$_D_x\n", obj.stage(1).bearing(1).C_x);
            fprintf(new_ID, "   subvar($_D_y, str = '%e N/(m/s)')                                                        ! $SVG_bearing_damping.$_D_y\n", obj.stage(1).bearing(1).C_y);
            fprintf(new_ID, "   subvar($_D_z, str = '%e N/(m/s)')                                                       ! $SVG_bearing_damping.$_D_z\n", obj.stage(1).bearing(1).C_z);
            fprintf(new_ID, "   subvar($_D_alpha, str = '0.0 Nm/(rad/s)')                                                   ! $SVG_bearing_damping.$_D_alpha\n");
            fprintf(new_ID, "   subvar($_D_beta, str = '%e Nm/(rad/s)')                                                   ! $SVG_bearing_damping.$_D_beta\n", obj.stage(1).bearing(1).C_beta);
            fprintf(new_ID, "   subvar($_D_gamma, str = '%e Nm/(rad/s)')                                                  ! $SVG_bearing_damping.$_D_gamma\n", obj.stage(1).bearing(1).C_gamma);
            fprintf(new_ID, "subvargroup.end (                   $SVG_bearing_damping          )                           ! $SVG_bearing_damping\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin (                 $SVG_mesh_stiffness           )                           ! $SVG_mesh_stiffness\n");
            fprintf(new_ID, "   subvar($_flag, str = '%i', discr.desc (   1) = 'Basic Gear Pair (FE204)', discr.desc (   2) = 'Gear Pair (FE225)', discr.str (   1) = '1', discr.str (   2) = '0') ! $SVG_mesh_stiffness.$_flag\n", mesh_flag);
            fprintf(new_ID, "   subvar($_linear, str = 'SWITCH($SVG_mesh_stiffness.$_flag)\n{\n   CASE 1: 0;\n   CASE 0: 1;\n   DEFAULT: 0;\n}') ! $SVG_mesh_stiffness.$_linear\n");
            fprintf(new_ID, "   subvar($_advanced, str = 'SWITCH($SVG_mesh_stiffness.$_flag)\n{\n   CASE 1: 1;\n   CASE 0: 0;\n   DEFAULT: 1;\n}') ! $SVG_mesh_stiffness.$_advanced\n");
            fprintf(new_ID, "subvargroup.end (                   $SVG_mesh_stiffness           )                           ! $SVG_mesh_stiffness\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "\n");
            
            fclose(new_ID);
            
        end
        
        function update_subvar_Nejad(obj, varargin)
            %UPDATE_SUBVAR update some values of the .Nejad.SubVar file 
            % from the SIMPACK model of the NREL_5MW Drivetrain.
            %
            % Kinematic Tree:
            % ================================================================================
            %
            % --  $R_Isys
            % --   ---->$B_Bed_Plate.................................... [6 dof]  Joint $J_Bed_Plate of type 20
            % .
            % $CTN_bed_plate_motion..................................... [6 cft]  Connection of type 61 (Constraint type 20)
            % OR (see $SVG_base_input.$_flag)
            % .
            % $CTN_bed_plate_fixed...................................... [6 cft]  Connection of type 2 (Constraint type 0)
            %           ---->$Gearbox_Frame............................. [0 dof]  Joint $J_Gearbox_Frame of type 0
            %                ---->$B_ISShaft............................ [6 dof]  Joint $J_ISShaft of type 20
            %                     ---->$BG_stage_01.$B_Sun.............. [0 dof]  Joint $BG_stage_01.$J_Sun of type 0
            %                     ---->$BG_stage_02.$B_PLC.............. [0 dof]  Joint $BG_stage_02.$J_PLC of type 0
            %                          ---->$BG_stage_02.$B_Pin1........ [6 dof]  Joint $BG_stage_02.$J_Pin1 of type 20
            %                               ---->$BG_stage_02.$B_PL1.... [0 dof]  Joint $BG_stage_02.$J_PL1 of type 0
            %
            %                          ---->$BG_stage_02.$B_Pin2........ [6 dof]  Joint $BG_stage_02.$J_Pin2 of type 20
            %                               ---->$BG_stage_02.$B_PL2.... [0 dof]  Joint $BG_stage_02.$J_PL2 of type 0
            %
            %                          ---->$BG_stage_02.$B_Pin3........ [6 dof]  Joint $BG_stage_02.$J_Pin3 of type 20
            %                               ---->$BG_stage_02.$B_PL3.... [0 dof]  Joint $BG_stage_02.$J_PL3 of type 0
            %
            %           ---->$B_main_shaft_LSS.......................... [6 dof]  Joint $J_main_shaft_LSS of type 20
            %                ---->$BG_stage_01.$B_PLC................... [0 dof]  Joint $BG_stage_01.$J_PLC of type 0
            %                     ---->$BG_stage_01.$B_Pin1............. [6 dof]  Joint $BG_stage_01.$J_Pin1 of type 20
            %                          ---->$BG_stage_01.$B_PL1......... [0 dof]  Joint $BG_stage_01.$J_PL1 of type 0
            %
            %                     ---->$BG_stage_01.$B_Pin2............. [6 dof]  Joint $BG_stage_01.$J_Pin2 of type 20
            %                          ---->$BG_stage_01.$B_PL2......... [0 dof]  Joint $BG_stage_01.$J_PL2 of type 0
            %
            %                     ---->$BG_stage_01.$B_Pin3............. [6 dof]  Joint $BG_stage_01.$J_Pin3 of type 20
            %                          ---->$BG_stage_01.$B_PL3......... [0 dof]  Joint $BG_stage_01.$J_PL3 of type 0
            %
            %                ---->$B_HSIShaft........................... [6 dof]  Joint $J_HSIShaft of type 20
            %                     ---->$BG_stage_02.$B_Sun.............. [0 dof]  Joint $BG_stage_02.$J_Sun of type 0
            %                     ---->$BG_stage_03.$B_Gear............. [0 dof]  Joint $BG_stage_03.$J_Gear of type 0
            %
            %                ---->$B_HSShaft............................ [6 dof]  Joint $J_HSShaft of type 20
            %                     ---->$BG_stage_03.$B_Pinion........... [0 dof]  Joint $BG_stage_03.$J_Pinion of type 0
            %
            %                ---->$BG_stage_01.$B_Ring.................. [0 dof]  Joint $BG_stage_01.$J_Ring of type 0
            %                ---->$BG_stage_02.$B_Ring.................. [0 dof]  Joint $BG_stage_02.$J_Ring of type 0
            %
            %
            % color for the different elements in SIMPACK:
            % sun / pinion  : 231-049-051
            % planet / wheel: 075-139-191
            % ring          : 095-183-092
            % carrier       : 255-140-025
            % shaft         : 221-207-110
            % bearing       : 175-103-061
            % highlight     : 248-142-197
            %
            
            default = {'linear_mesh', false, ...
                       'gen_mode'   , 1, ...
                       'bed_plate'  , true};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            gamma_x = mean([obj.gamma('m_n1') obj.gamma('m_n2') obj.gamma('m_n3')]);
            gamma_load = power(gamma_x, 2.0);
            gamma_Nm = power(gamma_x, 3.0);
            
%             k_sp = zeros(2, 1);
%             k_rp = k_sp;
%             for idx = 1:2
%                 k_sp(idx) = obj.stage(idx).sub_set('sun_planet').k_mesh;
%                 k_rp(idx) = obj.stage(idx).sub_set('planet_ring').k_mesh;
%             end
%             
%             k_3 = obj.stage(idx + 1).k_mesh;
            
            k_sp = [8.8807 5.6255]*1.0e9;
            k_rp = [1.006 2.88]*1.0e10;
            k_3 = 8.06e9;

            new_ID = fopen('@NREL_5MW\NREL_5MW.Nejad.subvar', 'w');
            
            fprintf(new_ID, "!file.version=3.3! Removing this line will make the file unreadable\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "!**********************************************************************\n");
            fprintf(new_ID, "! SubVars\n");
            fprintf(new_ID, "!**********************************************************************\n");
            fprintf(new_ID, "subvar($_gamma_Nm, str = '%g')                                        ! $_gamma_Nm\n", gamma_Nm);
            fprintf(new_ID, "subvar($_gamma_load, str = '%g')                                      ! $_gamma_load\n", gamma_load);
            fprintf(new_ID, "subvar($_rated_power, str = '%g kW')                                  ! $_rated_power\n", obj.P_rated);
            fprintf(new_ID, "subvar($_rotor_speed, str = '%g rpm')                                 ! $_rotor_speed\n", obj.n_rotor);
            fprintf(new_ID, "subvar($_GB_ratio, str = '%g')                                        ! $_GB_ratio\n", obj.u);
            fprintf(new_ID, "subvargroup.begin($SVG_main_shaft_LSS)                                ! $SVG_main_shaft_LSS\n");
            fprintf(new_ID, "   subvar($_gamma_d, str = '%g')                                      ! $SVG_main_shaft_LSS.$_gamma_d\n", obj.gamma('d_S'));
            fprintf(new_ID, "   subvar($_gamma_L, str = '%g')                                      ! $SVG_main_shaft_LSS.$_gamma_L\n", obj.gamma('L_S'));
            fprintf(new_ID, "   subvargroup.begin($SVG_INP_A)                                      ! $SVG_main_shaft_LSS.$SVG_INP_A\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '0 N/m')                                 ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '1.54e+10 N/m')                          ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '15400000000 N/m')                       ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                             ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '5e+6 Nm/rad')                          ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '5e+6 Nm/rad')                          ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '0 Ns/m')                                ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                        ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                        ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 N m s/rad')                          ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                    ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                    ! $SVG_main_shaft_LSS.$SVG_INP_A.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_INP_A)                                        ! $SVG_main_shaft_LSS.$SVG_INP_A\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_INP_B)                                      ! $SVG_main_shaft_LSS.$SVG_INP_B\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '4.06e+8 N/m')                           ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '1.54e+10 N/m')                          ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '15400000000 N/m')                       ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                             ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '97637.5 Nm/rad')                       ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '97637.5 Nm/rad')                       ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                        ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                        ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                        ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                            ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                    ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                    ! $SVG_main_shaft_LSS.$SVG_INP_B.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_INP_B)                                        ! $SVG_main_shaft_LSS.$SVG_INP_B\n");
            fprintf(new_ID, "subvargroup.end($SVG_main_shaft_LSS)                                  ! $SVG_main_shaft_LSS\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin($SVG_stage_01)  ! $SVG_stage_01\n");
            fprintf(new_ID, "   subvar($_gamma_m_n, str = '%g')                                    ! $SVG_stage_01.$_gamma_m_n\n", obj.gamma('m_n1'));
            fprintf(new_ID, "   subvar($_gamma_b, str = '%g')                                      ! $SVG_stage_01.$_gamma_b\n", obj.gamma('b_1'));
            fprintf(new_ID, "   subvar($_gamma_d, str = '%g')                                      ! $SVG_stage_01.$_gamma_d\n", obj.gamma('d_1'));
            fprintf(new_ID, "   subvar($_gamma_L, str = '%g')                                      ! $SVG_stage_01.$_gamma_L\n", obj.gamma('L_1'));
            fprintf(new_ID, "   subvar($_m_n, str = '%g mm')                                       ! $SVG_stage_01.$_m_n\n", obj.stage(1).m_n);
            fprintf(new_ID, "   subvar($_b, str = '%g mm')                                         ! $SVG_stage_01.$_b\n", obj.stage(1).b);
            fprintf(new_ID, "   subvar($_a_w, str = '%g mm')                                       ! $SVG_stage_01.$_a_w\n", obj.stage(1).a_w);
            fprintf(new_ID, "   subvar($_alpha, str = '20 deg')                                    ! $SVG_stage_01.$_alpha\n");
            fprintf(new_ID, "   subvar($_x_sun, str = '0.617')                                     ! $SVG_stage_01.$_x_sun\n");
            fprintf(new_ID, "   subvar($_x_planet, str = '0.802')                                  ! $SVG_stage_01.$_x_planet\n");
            fprintf(new_ID, "   subvar($_x_ring, str = '-0.501')                                   ! $SVG_stage_01.$_x_ring\n");
            fprintf(new_ID, "   subvar($_z_sun, str = '19')                                        ! $SVG_stage_01.$_z_sun\n");
            fprintf(new_ID, "   subvar($_z_planet, str = '17')                                     ! $SVG_stage_01.$_z_planet\n");
            fprintf(new_ID, "   subvar($_z_ring, str = '56')                                       ! $SVG_stage_01.$_z_ring\n");
            fprintf(new_ID, "   subvar($_Norm_Backlash, str = '0.40*$SVG_stage_01.$_gamma_m_n mm') ! $SVG_stage_01.$_Norm_Backlash\n");
            fprintf(new_ID, "   subvar($_Tooth_Stiff_ratio, str = '0.8')                           ! $SVG_stage_01.$_Tooth_Stiff_ratio\n");
            fprintf(new_ID, "   subvar($_k_SP, str = '%g N/m')                                     ! $SVG_stage_01.$_k_SP\n", k_sp(1));
            fprintf(new_ID, "   subvar($_k_RP, str = '%g N/m')                                     ! $SVG_stage_01.$_k_RP\n", k_rp(1));
            fprintf(new_ID, "   subvar($_mesh_damping, str = '5.0e8 Ns/m')                         ! $SVG_stage_01.$_mesh_damping\n");
            fprintf(new_ID, "   subvar($_transition_depth, str = '$SVG_mesh.$_full_damping')       ! $SVG_stage_01.$_transition_depth\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_PLC_A)                                      ! $SVG_stage_01.$SVG_PLC_A\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '6.6e+4 N/m')                            ! $SVG_stage_01.$SVG_PLC_A.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '1770000000 N/m')                        ! $SVG_stage_01.$SVG_PLC_A.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '1770000000 N/m')                        ! $SVG_stage_01.$SVG_PLC_A.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                             ! $SVG_stage_01.$SVG_PLC_A.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '1390000 Nm/rad')                       ! $SVG_stage_01.$SVG_PLC_A.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '1.39e+6 Nm/rad')                       ! $SVG_stage_01.$SVG_PLC_A.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                        ! $SVG_stage_01.$SVG_PLC_A.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                        ! $SVG_stage_01.$SVG_PLC_A.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                        ! $SVG_stage_01.$SVG_PLC_A.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                            ! $SVG_stage_01.$SVG_PLC_A.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                    ! $SVG_stage_01.$SVG_PLC_A.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                    ! $SVG_stage_01.$SVG_PLC_A.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_PLC_A)                                        ! $SVG_stage_01.$SVG_PLC_A\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_PLC_B)  ! $SVG_stage_01.$SVG_PLC_B\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '6.6e+7 N/m')                                ! $SVG_stage_01.$SVG_PLC_B.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '1770000000 N/m')                            ! $SVG_stage_01.$SVG_PLC_B.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '1770000000 N/m')                            ! $SVG_stage_01.$SVG_PLC_B.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_01.$SVG_PLC_B.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '1390000 Nm/rad')                              ! $SVG_stage_01.$SVG_PLC_B.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '1.39e+6 Nm/rad')                              ! $SVG_stage_01.$SVG_PLC_B.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_01.$SVG_PLC_B.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_01.$SVG_PLC_B.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_01.$SVG_PLC_B.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_01.$SVG_PLC_B.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_01.$SVG_PLC_B.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_01.$SVG_PLC_B.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_PLC_B)  ! $SVG_stage_01.$SVG_PLC_B\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_PL_AB)  ! $SVG_stage_01.$SVG_PL_AB\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '91016 N/m')                                 ! $SVG_stage_01.$SVG_PL_AB.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '9.37E+9 N/m')                               ! $SVG_stage_01.$SVG_PL_AB.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '3.20E+9 N/m')                               ! $SVG_stage_01.$SVG_PL_AB.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_01.$SVG_PL_AB.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '1.39E+6 Nm/rad')                              ! $SVG_stage_01.$SVG_PL_AB.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '4.45E+6 Nm/rad')                              ! $SVG_stage_01.$SVG_PL_AB.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_01.$SVG_PL_AB.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_01.$SVG_PL_AB.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_01.$SVG_PL_AB.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_01.$SVG_PL_AB.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_01.$SVG_PL_AB.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_01.$SVG_PL_AB.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_PL_AB)  ! $SVG_stage_01.$SVG_PL_AB\n");
            fprintf(new_ID, "subvargroup.end($SVG_stage_01)  ! $SVG_stage_01\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin($SVG_stage_02)  ! $SVG_stage_02\n");
            fprintf(new_ID, "   subvar($_gamma_m_n, str = '%g')                                     ! $SVG_stage_02.$_gamma_m_n\n", obj.gamma('m_n2'));
            fprintf(new_ID, "   subvar($_gamma_b, str = '%g')                                       ! $SVG_stage_02.$_gamma_b\n", obj.gamma('b_2'));
            fprintf(new_ID, "   subvar($_gamma_d, str = '%g')                                       ! $SVG_stage_02.$_gamma_d\n", obj.gamma('d_2'));
            fprintf(new_ID, "   subvar($_gamma_L, str = '%g')                                       ! $SVG_stage_02.$_gamma_L\n", obj.gamma('L_2'));
            fprintf(new_ID, "   subvar($_m_n, str = '%g mm')                                       ! $SVG_stage_02.$_m_n\n", obj.stage(2).m_n);
            fprintf(new_ID, "   subvar($_b, str = '%g mm')                                        ! $SVG_stage_02.$_b\n", obj.stage(2).b);
            fprintf(new_ID, "   subvar($_a_w, str = '%g mm')                                      ! $SVG_stage_02.$_a_w\n", obj.stage(2).a_w);
            fprintf(new_ID, "   subvar($_alpha, str = '20 deg')                                    ! $SVG_stage_02.$_alpha\n");
            fprintf(new_ID, "   subvar($_x_sun, str = '0.389')                                     ! $SVG_stage_02.$_x_sun\n");
            fprintf(new_ID, "   subvar($_x_planet, str = '0.504')                                  ! $SVG_stage_02.$_x_planet\n");
            fprintf(new_ID, "   subvar($_x_ring, str = '0.117')                                    ! $SVG_stage_02.$_x_ring\n");
            fprintf(new_ID, "   subvar($_z_sun, str = '18')                                        ! $SVG_stage_02.$_z_sun\n");
            fprintf(new_ID, "   subvar($_z_planet, str = '36')                                     ! $SVG_stage_02.$_z_planet\n");
            fprintf(new_ID, "   subvar($_z_ring, str = '93')                                       ! $SVG_stage_02.$_z_ring\n");
            fprintf(new_ID, "   subvar($_Norm_Backlash, str = '0.40*$SVG_stage_02.$_gamma_m_n mm') ! $SVG_stage_02.$_Norm_Backlash\n");
            fprintf(new_ID, "   subvar($_Tooth_Stiff_ratio, str = '0.8')                           ! $SVG_stage_02.$_Tooth_Stiff_ratio\n");
            fprintf(new_ID, "   subvar($_k_SP, str = '%g N/m')                                  ! $SVG_stage_02.$_k_SP\n", k_sp(2));
            fprintf(new_ID, "   subvar($_k_RP, str = '%g N/m')                                  ! $SVG_stage_02.$_k_RP\n", k_rp(2));
            fprintf(new_ID, "   subvar($_mesh_damping, str = '5.0e8 Ns/m')                         ! $SVG_stage_02.$_mesh_damping\n");
            fprintf(new_ID, "   subvar($_transition_depth, str = '$SVG_mesh.$_full_damping')       ! $SVG_stage_02.$_transition_depth\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_IMS_PLC_A)  ! $SVG_stage_02.$SVG_IMS_PLC_A\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '9.1e+4 N/m')                                ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '6.12e+7 N/m')                               ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '1.16e+9 N/m')                               ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '74428 Nm/rad')                                ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '3239.8 Nm/rad')                               ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_02.$SVG_IMS_PLC_A.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_IMS_PLC_A)  ! $SVG_stage_02.$SVG_IMS_PLC_A\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_IMS_PLC_B)  ! $SVG_stage_02.$SVG_IMS_PLC_B\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '9.1e+7 N/m')                                ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '6.12e+7 N/m')                               ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '1.16e+9 N/m')                               ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '74428 Nm/rad')                                ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '3239.8 Nm/rad')                               ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_02.$SVG_IMS_PLC_B.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_IMS_PLC_B)  ! $SVG_stage_02.$SVG_IMS_PLC_B\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_PL_AB)  ! $SVG_stage_02.$SVG_PL_AB\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '91016 N/m')                                 ! $SVG_stage_02.$SVG_PL_AB.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '6.12E+7 N/m')                               ! $SVG_stage_02.$SVG_PL_AB.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '1.16E+9 N/m')                               ! $SVG_stage_02.$SVG_PL_AB.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_02.$SVG_PL_AB.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '7.5e+4 Nm/rad')                               ! $SVG_stage_02.$SVG_PL_AB.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '7.5e+4 Nm/rad')                               ! $SVG_stage_02.$SVG_PL_AB.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_02.$SVG_PL_AB.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_02.$SVG_PL_AB.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_02.$SVG_PL_AB.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_02.$SVG_PL_AB.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_02.$SVG_PL_AB.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_02.$SVG_PL_AB.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_PL_AB)  ! $SVG_stage_02.$SVG_PL_AB\n");
            fprintf(new_ID, "subvargroup.end($SVG_stage_02)  ! $SVG_stage_02\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin($SVG_stage_03)  ! $SVG_stage_03\n");
            fprintf(new_ID, "   subvar($_gamma_m_n, str = '%g')                                     ! $SVG_stage_03.$_gamma_m_n\n", obj.gamma('m_n3'));
            fprintf(new_ID, "   subvar($_gamma_b, str = '%g')                                       ! $SVG_stage_03.$_gamma_b\n", obj.gamma('b_3'));
            fprintf(new_ID, "   subvar($_gamma_d, str = '%g')                                       ! $SVG_stage_03.$_gamma_d\n", obj.gamma('d_3'));
            fprintf(new_ID, "   subvar($_gamma_L, str = '%g')                                       ! $SVG_stage_03.$_gamma_L\n", obj.gamma('L_3'));
            fprintf(new_ID, "   subvar($_m_n, str = '%g mm')                                       ! $SVG_stage_03.$_m_n\n", obj.stage(3).m_n);
            fprintf(new_ID, "   subvar($_b, str = '%g mm')                                        ! $SVG_stage_03.$_b\n", obj.stage(3).b);
            fprintf(new_ID, "   subvar($_a_w, str = '%g mm')                                      ! $SVG_stage_03.$_a_w\n", obj.stage(3).a_w);
            fprintf(new_ID, "   subvar($_alpha, str = '20 deg')                                    ! $SVG_stage_03.$_alpha\n");
            fprintf(new_ID, "   subvar($_helix, str = '10 deg')                                    ! $SVG_stage_03.$_helix\n");
            fprintf(new_ID, "   subvar($_x_pinion, str = '0.480')                                  ! $SVG_stage_03.$_x_pinion\n");
            fprintf(new_ID, "   subvar($_x_gear, str = '0.669')                                    ! $SVG_stage_03.$_x_gear\n");
            fprintf(new_ID, "   subvar($_z_pinion, str = '24')                                     ! $SVG_stage_03.$_z_pinion\n");
            fprintf(new_ID, "   subvar($_z_gear, str = '95')                                       ! $SVG_stage_03.$_z_gear\n");
            fprintf(new_ID, "   subvar($_Norm_Backlash, str = '0.40*$SVG_stage_03.$_gamma_m_n mm') ! $SVG_stage_03.$_Norm_Backlash\n");
            fprintf(new_ID, "   subvar($_Tooth_Stiff_ratio, str = '0.8')                           ! $SVG_stage_03.$_Tooth_Stiff_ratio\n");
            fprintf(new_ID, "   subvar($_k_mesh, str = '%g N/m')                                ! $SVG_stage_03.$_k_mesh\n", k_3);
            fprintf(new_ID, "   subvar($_mesh_damping, str = '5.0e8 Ns/m')                         ! $SVG_stage_03.$_mesh_damping\n");
            fprintf(new_ID, "   subvar($_transition_depth, str = '$SVG_mesh.$_full_damping')       ! $SVG_stage_03.$_transition_depth\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_IMS_A)  ! $SVG_stage_03.$SVG_IMS_A\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '0 N/m')                                     ! $SVG_stage_03.$SVG_IMS_A.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '6.12E+7 N/m')                               ! $SVG_stage_03.$SVG_IMS_A.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '1.16E+9 N/m')                               ! $SVG_stage_03.$SVG_IMS_A.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_03.$SVG_IMS_A.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '74428 Nm/rad')                                ! $SVG_stage_03.$SVG_IMS_A.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '3239.8 Nm/rad')                               ! $SVG_stage_03.$SVG_IMS_A.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_A.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_A.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_A.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_03.$SVG_IMS_A.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_IMS_A.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_IMS_A.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_IMS_A)  ! $SVG_stage_03.$SVG_IMS_A\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_IMS_B)  ! $SVG_stage_03.$SVG_IMS_B\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '7.41e+7 N/m')                               ! $SVG_stage_03.$SVG_IMS_B.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '5.17E+8 N/m')                               ! $SVG_stage_03.$SVG_IMS_B.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '484000000 N/m')                             ! $SVG_stage_03.$SVG_IMS_B.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_03.$SVG_IMS_B.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '1.67e+6 Nm/rad')                              ! $SVG_stage_03.$SVG_IMS_B.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '1.78e+6 Nm/rad')                              ! $SVG_stage_03.$SVG_IMS_B.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_B.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_B.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_B.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_03.$SVG_IMS_B.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_IMS_B.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_IMS_B.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_IMS_B)  ! $SVG_stage_03.$SVG_IMS_B\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_IMS_C)  ! $SVG_stage_03.$SVG_IMS_C\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '78723000 N/m')                              ! $SVG_stage_03.$SVG_IMS_C.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '7.3708E+8 N/m')                             ! $SVG_stage_03.$SVG_IMS_C.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '3.2612e+8 N/m')                             ! $SVG_stage_03.$SVG_IMS_C.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_03.$SVG_IMS_C.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '1123400 Nm/rad')                              ! $SVG_stage_03.$SVG_IMS_C.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '2539000 Nm/rad')                              ! $SVG_stage_03.$SVG_IMS_C.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_C.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_C.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_03.$SVG_IMS_C.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_03.$SVG_IMS_C.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_IMS_C.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_IMS_C.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_IMS_C)  ! $SVG_stage_03.$SVG_IMS_C\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_HS_A )  ! $SVG_stage_03.$SVG_HS_A\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '1.26e+8 N/m')                               ! $SVG_stage_03.$SVG_HS_A.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '8.21e+8 N/m')                               ! $SVG_stage_03.$SVG_HS_A.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '821000000 N/m')                             ! $SVG_stage_03.$SVG_HS_A.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_03.$SVG_HS_A.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '1123400 Nm/rad')                              ! $SVG_stage_03.$SVG_HS_A.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '2539000 Nm/rad')                              ! $SVG_stage_03.$SVG_HS_A.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_03.$SVG_HS_A.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_03.$SVG_HS_A.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_03.$SVG_HS_A.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_03.$SVG_HS_A.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_HS_A.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_HS_A.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_HS_A )  ! $SVG_stage_03.$SVG_HS_A\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_HS_B )  ! $SVG_stage_03.$SVG_HS_B\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '6.7e+7 N/m')                                ! $SVG_stage_03.$SVG_HS_B.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '8.09e+8 N/m')                               ! $SVG_stage_03.$SVG_HS_B.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '1.33e+8 N/m')                               ! $SVG_stage_03.$SVG_HS_B.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_03.$SVG_HS_B.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '1.71e+5 Nm/rad')                              ! $SVG_stage_03.$SVG_HS_B.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '1.04e+6 Nm/rad')                              ! $SVG_stage_03.$SVG_HS_B.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_03.$SVG_HS_B.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_03.$SVG_HS_B.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_03.$SVG_HS_B.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_03.$SVG_HS_B.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_HS_B.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_HS_B.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_HS_B)  ! $SVG_stage_03.$SVG_HS_B\n");
            fprintf(new_ID, "   subvargroup.begin($SVG_HS_C)  ! $SVG_stage_03.$SVG_HS_C\n");
            fprintf(new_ID, "      subvar($_PL_c_x, str = '7.93e+7 N/m')                               ! $SVG_stage_03.$SVG_HS_C.$_PL_c_x\n");
            fprintf(new_ID, "      subvar($_PL_c_y, str = '1.04e+9 N/m')                               ! $SVG_stage_03.$SVG_HS_C.$_PL_c_y\n");
            fprintf(new_ID, "      subvar($_PL_c_z, str = '7.29e+7 N/m')                               ! $SVG_stage_03.$SVG_HS_C.$_PL_c_z\n");
            fprintf(new_ID, "      subvar($_PL_c_al, str = '0 Nm/rad')                                    ! $SVG_stage_03.$SVG_HS_C.$_PL_c_al\n");
            fprintf(new_ID, "      subvar($_PL_c_be, str = '97637.5 Nm/rad')                              ! $SVG_stage_03.$SVG_HS_C.$_PL_c_be\n");
            fprintf(new_ID, "      subvar($_PL_c_ga, str = '1.39e+6 Nm/rad')                              ! $SVG_stage_03.$SVG_HS_C.$_PL_c_ga\n");
            fprintf(new_ID, "      subvar($_PL_d_x, str = '4.53e+002 Ns/m')                             ! $SVG_stage_03.$SVG_HS_C.$_PL_d_x\n");
            fprintf(new_ID, "      subvar($_PL_d_y, str = '4.20e+004 Ns/m')                             ! $SVG_stage_03.$SVG_HS_C.$_PL_d_y\n");
            fprintf(new_ID, "      subvar($_PL_d_z, str = '3.06e+004 Ns/m')                             ! $SVG_stage_03.$SVG_HS_C.$_PL_d_z\n");
            fprintf(new_ID, "      subvar($_PL_d_al, str = '0 Nms/rad')                                    ! $SVG_stage_03.$SVG_HS_C.$_PL_d_al\n");
            fprintf(new_ID, "      subvar($_PL_d_be, str = '3.43e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_HS_C.$_PL_d_be\n");
            fprintf(new_ID, "      subvar($_PL_d_ga, str = '4.78e+001 Nms/rad')                            ! $SVG_stage_03.$SVG_HS_C.$_PL_d_ga\n");
            fprintf(new_ID, "   subvargroup.end($SVG_HS_C)  ! $SVG_stage_03.$SVG_HS_C\n");
            fprintf(new_ID, "subvargroup.end($SVG_stage_03)  ! $SVG_stage_03\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin($SVG_mesh)                           ! $SVG_mesh\n");
            fprintf(new_ID, "   subvar($_full_damping, str= '10e-6 m')                                                     ! $SVG_mesh.$_full_damping\n");
            fprintf(new_ID, "   subvar($_flag, str = '%i', discr.desc(1) = 'detailed (FE225)', discr.desc(2) = 'linear (FE204)', discr.str (   1) = '0', discr.str (   2) = '1') ! $SVG_mesh.$_flag\n", default.linear_mesh);
            fprintf(new_ID, "   subvar($_use_detailed, str= 'IF($SVG_mesh.$_flag == 1)\n{\n\t1\n}\nELSE\n{\n\t0\n}')                      ! $SVG_mesh.$_use_detailed\n");
            fprintf(new_ID, "   subvar($_use_linear,   str= 'IF($SVG_mesh.$_flag == 1)\n{\n\t0\n}\nELSE\n{\n\t1\n}')                      ! $SVG_mesh.$_use_linear\n");
            fprintf(new_ID, "subvargroup.end($SVG_mesh)                           ! $SVG_mesh\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin($SVG_generator)                           ! $SVG_generator\n");
            fprintf(new_ID, "   subvar($_flag, str = '%i', discr.desc(1) = 'torque from speed, T(omega)', discr.desc(2) = 'PI control (CE129)', discr.desc(3) = 'OFF', discr.str (   1) = '1', discr.str (   2) = '2', discr.str(3) = '0') ! $SVG_generator.$_flag\n", default.gen_mode);
            fprintf(new_ID, "   subvar($_use_torque,  str = 'IF($SVG_generator.$_flag == 1)\n{\n\t0\n}\nELSE\n{\n\t1\n}') ! $SVG_generator.$_use_torque\n");
            fprintf(new_ID, "   subvar($_use_control, str = 'IF($SVG_generator.$_flag == 2)\n{\n\t0\n}\nELSE\n{\n\t1\n}') ! $SVG_generator.$_use_control\n");
            fprintf(new_ID, "subvargroup.end($SVG_generator)                           ! $SVG_generator\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "subvargroup.begin($SVG_base_input)                           ! $SVG_base_input\n");
            fprintf(new_ID, "   subvar($_flag, str = '%i', discr.desc(1) = 'ON', discr.desc(2) = 'OFF', discr.str(1) = '1', discr.str(2) = '0') ! $SVG_base_input.$_flag\n", default.bed_plate);
            fprintf(new_ID, "   subvar($_ON,  str= 'IF($SVG_base_input.$_flag == 1)\n{\n\t1\n}\nELSE\n{\n\t0\n}')                       ! $SVG_base_input.$_ON\n");
            fprintf(new_ID, "   subvar($_OFF, str= 'IF($SVG_base_input.$_flag == 1)\n{\n\t0\n}\nELSE\n{\n\t1\n}')                       ! $SVG_base_input.$_OFF\n");
            fprintf(new_ID, "subvargroup.end($SVG_base_input)                           ! $SVG_base_input\n");
            fprintf(new_ID, "\n");
            fprintf(new_ID, "\n");
            
            fclose(new_ID);
            
        end
        
        function update_subvar_WTDB(obj, varargin)
            %UPDATE_SUBVAR_WTDB updates some parameters from Simpack's Wind
            % Turbine database gearbox.
            % 
            % Kinematic Tree:
            % ================================================================================
            %   $R_Isys
            %    ---->$B_HSG................... [6 dof]  Connection $G_SHAFT.$CTN_HSG of type 63 (Joint type 15)
            %         ---->$B_LSS_PLC.......... [6 dof]  Connection $G_SHAFT.$CTN_LSS_PLC of type 63 (Joint type 15)
            %              ---->$B_LSS_PLT1.... [6 dof]  Connection $G_SHAFT.$CTN_LSS_PLT1 of type 63 (Joint type 15)
            %              ---->$B_LSS_PLT2.... [6 dof]  Connection $G_SHAFT.$CTN_LSS_PLT2 of type 63 (Joint type 15)
            %              ---->$B_LSS_PLT3.... [6 dof]  Connection $G_SHAFT.$CTN_LSS_PLT3 of type 63 (Joint type 15)
            %              ---->$B_LSS_PLT4.... [6 dof]  Connection $G_SHAFT.$CTN_LSS_PLT4 of type 63 (Joint type 15)
            % 
            %         ---->$B_IMS_PLC.......... [6 dof]  Connection $G_SHAFT.$CTN_IMS_PLC of type 63 (Joint type 15)
            %              ---->$B_LSS_SUN..... [6 dof]  Connection $G_SHAFT.$CTN_LSS_SUN of type 63 (Joint type 15)
            %              ---->$B_IMS_PLT1.... [6 dof]  Connection $G_SHAFT.$CTN_IMS_PLT1 of type 63 (Joint type 15)
            %              ---->$B_IMS_PLT2.... [6 dof]  Connection $G_SHAFT.$CTN_IMS_PLT2 of type 63 (Joint type 15)
            %              ---->$B_IMS_PLT3.... [6 dof]  Connection $G_SHAFT.$CTN_IMS_PLT3 of type 63 (Joint type 15)
            % 
            %         ---->$B_HSS_BGR.......... [6 dof]  Connection $G_SHAFT.$CTN_HSS_BGR of type 63 (Joint type 15)
            %              ---->$B_IMS_SUN..... [6 dof]  Connection $G_SHAFT.$CTN_IMS_SUN of type 63 (Joint type 15)
            % 
            %         ---->$B_HSS_HSH.......... [6 dof]  Connection $G_SHAFT.$CTN_HSS_HSH of type 63 (Joint type 15)
            %
            
        end
        
    end
    
end
