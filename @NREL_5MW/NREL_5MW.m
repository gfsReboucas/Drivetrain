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
            gamma = {'P'   , 1.0, 'n_R' , 1.0, ...
                     'm_n1', 1.0, 'b_1' , 1.0, 'd_1' , 1.0, 'L_1' , 1.0, ...
                     'm_n2', 1.0, 'b_2' , 1.0, 'd_2' , 1.0, 'L_2' , 1.0, ...
                     'm_n3', 1.0, 'b_3' , 1.0, 'd_3' , 1.0, 'L_3' , 1.0, ...
                     'M_R' , 1.0, 'J_R' , 1.0, ...
                     'M_G' , 1.0, 'J_G' , 1.0, ...
                     'd_s' , 1.0, 'L_s' , 1.0, ...
                     'd_S' , 1.0, 'L_S' , 1.0, ...
                     'fault_type'   ,  '', ...
                     'fault_val'    , 0.0, ...
                     'dynamic_model', @Kahraman_94};
%                      'dynamic_model', @Lin_Parker_99};
%                      'dynamic_model', @Dynamic_Formulation};
            
            gamma = scaling_factor.process_varargin(gamma, varargin);
            
%             dynamic_model = gamma.dynamic_model; % No fault
            dynamic_model = gamma.dynamic_model;
            fault_type    = gamma.fault_type;
            fault_val     = gamma.fault_val;
            
            gamma = rmfield(gamma, {'dynamic_model', 'fault_type', 'fault_val'});
            gamma = scaling_factor(gamma);

            N_st = 3;
            stage = repmat(Gear_Set(), 1, N_st);
            
            for idx = 1:N_st
                stg = NREL_5MW.gear_stage(idx);
                
                gamma_idx = gamma.ends_with(num2str(idx));
                stage(idx) = stg.scale_by(gamma_idx);
            end
            
            P_r = gamma('P')* 5.0e3; % [kW]    , Rated power
            n_r = gamma('n_R')*12.1; % [1/min.], Input speed
            
            LSS = NREL_5MW.shaft(0);
            
            inp_shaft = Shaft('d'       , LSS.d*gamma('d_s'), ...
                              'L'       , LSS.L*gamma('L_s'), ...
                              'bearing' , LSS.bearing, ...
                              'material', LSS.material);
            
            m_R =      110.0e3*gamma('M_R');
            J_R = 57231535.0  *gamma('J_R'); % according to [1, 3]
            m_G =     1900.0  *gamma('M_G');
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
                           'dynamic_model', dynamic_model, ...
                           'fault_type'   , fault_type, ...
                           'fault_val'    , fault_val);

            obj.gamma = scaling_factor(gamma.name, gamma.value);

        end
    end
    
    %% Calculation:
    methods
        function [SH, SF, SShaft, calc] = safety_factor_stage(obj, idx)
            L_h    =  20.0  * ... % year
                     365.25 * ... % day/year
                      24.0;       % hours/day
            K_A    = 1.25;
            
            switch(idx)
                case 1
                    n_idx = [nan, nan, 0.0, obj.n_rotor];
                    extension = 'Z14';
                case 2
                    n_idx = [nan, nan, 0.0, obj.n_out(1)];
                    extension = 'Z14';
                case 3
                    n_idx = [obj.n_out(3), nan];
                    extension = 'Z12';
                otherwise
                    error('NREL_5MW:safety_factor', 'Index out of range.');
            end
            
            KS_file = sprintf('\\@%s\\stage_%02d.%s', class(obj), idx, extension);
            
%             calc = KSISO_6336(obj.stage(idx),  'P_rated'     , obj.P_rated, ...
            calc = MATISO_6336(obj.stage(idx), 'P_rated'     , obj.P_rated, ...
                                               'n_out'       , obj.n_out(idx), ...
                                               'S_Hmin'      , obj.S_Hmin, ...
                                               'S_Fmin'      , obj.S_Fmin, ...
                                               'L_h'         , L_h, ... % [h], Required life
                                               'K_A'         , K_A, ... % [-], Application factor
                                               'lubricant_ID', 11170, ...
                                               'nu_40'       , 220.0, ... % [], lubricant viscosity
                                               'n_nominal'   , n_idx, ...
                                               'KS_file'     , KS_file);
            [SH, SF] = calc.safety_factors('nu_40'       , 220.0, ...
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
            % - https://www.skf.com/group/products/rolling-bearings
            % - http://bearingsize.info/
            % - https://simplybearings.co.uk/shop/
            %
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    Name   |           Label          | a [-] | C [kN] | e [-] |  X [-] |  Y [-]  | K_x (N/m) | K_y (N/m) | K_z (N/m) | K_beta (Nm/rad) | K_gamma (Nm/rad) | Type |   D  |   d  |  B  |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   INP-A   |     SKF C 30/1250 MB*    |  10/3 |  22000 |   -   |    1   |    0    |    0.0    |   1.5e10  |   1.5e10  |      5.0e6      |       5.0e6      | CARB | 1750 | 1250 | 375 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   INP-B   |    SKF 231/750 CA/W33    |  10/3 |  16518 |  0.28 | 1/0.67 | 2.4/3.6 |   4.06e8  |  1.54e10  |  1.54e10  |       0.0       |        0.0       |  SRB | 1220 |  750 | 365 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   PLC-A   |    SKF 240/710 ECA/W33   |  10/3 |  11164 |  0.27 | 1/0.67 | 2.5/3.7 |   6.6e4   |   1.7e9   |   1.1e9   |      5.6e5      |       1.3e5      |  SRB | 1030 |  710 | 315 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   PLC-B   | SKF NF 28/1000 ECMP/HA1+ |  10/3 |  3690  |  0.3  | 1/0.92 |  0/0.4  |   6.6e7   |   1.7e9   |   1.1e9   |      5.6e5      |       1.3e5      |  CRB | 1220 | 1000 | 128 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    PL-A   |     SKF NNCF 5080 CV     |  10/3 |  5500  |  0.15 | 1/0.92 |  0/0.4  |   9.1e4   |   9.4e9   |   3.2e9   |      1.4e6      |       4.5e6      |  CRB |  600 |  400 | 272 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    PL-B   |     SKF NNCF 5080 CV     |  10/3 |  5500  |  0.15 | 1/0.92 |  0/0.4  |   9.1e4   |   9.4e9   |   3.2e9   |      1.4e6      |       4.5e6      |  CRB |  600 |  400 | 272 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % | IMS-PLC-A |       SKF C30/710M       |  10/3 |  8800  |   -   |    1   |    0    |   9.1e4   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      | CARB | 1030 |  710 | 236 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % | IMS-PLC-B |    SKF NU 20/600 ECMA    |  10/3 |  4180  |  0.3  | 1/0.92 |  0/0.4  |   9.1e7   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  870 |  600 | 155 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |  IMS-PL-A |     SKF NNCF 4976 CV+    |  10/3 |  2380  |  0.15 | 1/0.92 |  0/0.4  |   9.1e4   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  520 |  380 | 140 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |  IMS-PL-B |     SKF NNCF 4976 CV+    |  10/3 |  2380  |  0.15 | 1/0.92 |  0/0.4  |   9.1e4   |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  520 |  380 | 140 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   IMS-A   |     SKF NNCF 4880 CV+    |  10/3 |  1420  |  0.15 | 1/0.92 |  0/0.4  |     0     |    6e7    |   1.2e9   |      7.5e4      |       7.5e4      |  CRB |  500 |  400 | 100 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   IMS-B   |  SKF LM 567949/910/HA1 & |  10/3 |  1467  |  0.4  |  1/0.4 |  0/1.5  |   7.8e7   |   7.4e8   |   3.3e8   |      1.1e6      |       2.5e6      |  TRB |  550 |  410 |  86 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |   IMS-C   |  SKF LM 567949/910/HA1 & |  10/3 |  1467  |  0.4  |  1/0.4 |  0/1.5  |   7.8e7   |   7.4e8   |   3.3e8   |      1.1e6      |       2.5e6      |  TRB |  550 |  410 |  86 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    HS-A   |     SKF NCF 2240 ECJB    |  10/3 |  1460  |  0.3  | 1/0.92 |  0/0.4  |   1.3e8   |   8.2e8   |   8.2e8   |      1.7e5      |        1e6       |  CRB |  360 |  200 |  98 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    HS-B   |        SKF 32060 X ?     |  10/3 |  1644  |  0.43 |  1/0.4 |  0/1.4  |   6.7e7   |    8e8    |   1.3e8   |      1.7e5      |        1e6       |  TRB |  460 | 200? | 100 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % |    HS-C   |        SKF 32060 X ?     |  10/3 |  1644  |  0.43 |  1/0.4 |  0/1.4  |    8e7    |    1e9    |   7.3e7   |      1.7e5      |        1e6       |  TRB |  460 | 200? | 100 |
            % +-----------+--------------------------+-------+--------+-------+--------+---------+-----------+-----------+-----------+-----------------+------------------+------+------+------+-----+
            % *: not found on SKF's website (access 12.11.2020)
            % +: multiple bearings with the same geometric parameters
            % &: in imperial units and rounded values in mm
            % ?: suspect of a typo, the correct would be 300 instead of 200
            %
            
            cx = 453.0;         cy = 42.0e3;            cz = 30600.0;
                                cb = 34.3;              cg = 47.8;
            
            a = 10.0/3.0;       L_10 =  20.0  * ... % year
                                       365.25 * ... % day/year
                                        24.0;       % hours/day
                  
            switch(idx)
                case 0 % Main shaft
                    INP_A     = Bearing('name', 'INP-A', 'type'  ,  'CARB', 'x'      , 333.0, ...
                                        'a'   ,       a, 'C'     ,    2200.0e3, 'e'      , +inf, ...
                                        'X'   , [1.0 0.0], 'Y'     , [0.0 0.0], 'L_10'   , L_10, ...
                                        'K_x' ,   0.0e4, 'K_y'   ,  1.5e10, 'K_z'    , 1.5e10, ...
                                                         'K_beta',   5.0e6, 'K_gamma',  5.0e6, ...
                                        'C_x' ,  0.0*cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'D'   ,  1750.0, 'd'     ,  1250.0, 'B'      ,  375.0);
                    INP_B     = Bearing('name', 'INP-B', 'type'  , 'SRB'  , 'x'      , 1166.0, ...
                                        'a'   ,       a, 'C'     , 16518.0e3, 'e'      , 0.28, ...
                                        'X'   , [1.0 0.67], 'Y'  , [2.4 3.6], 'L_10'   , L_10, ...
                                        'K_x' ,  4.06e8, 'K_y'   , 1.54e10, 'K_z'    , 1.54e10, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'D'   ,  1220.0, 'd'     ,   750.0, 'B'      ,   365.0);
                    
                    brg = [INP_A, INP_B];
                    
                case 1
                    PL_A      = Bearing('name',  'PL-A', 'type'  ,  'CRB', ...
                                        'a'   ,       a, 'C'     , 5500.0e3, 'e'      , 0.15, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x' ,   9.1e4, 'K_y'   ,  9.4e9, 'K_z'    , 3.2e9, ...
                                                         'K_beta',  1.4e6, 'K_gamma', 4.5e6, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'D'   ,   600.0, 'd'     ,  400.0, 'B'      , 272.0);
                    PL_B      = Bearing('name',  'PL-B', 'type'  ,  'CRB', ...
                                        'a'   ,       a, 'C'     , 5500.0e3, 'e'      , 0.15, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x',    9.1e4, 'K_y'   ,  9.4e9, 'K_z'    , 3.2e9, ...
                                                         'K_beta',  1.4e6, 'K_gamma', 4.5e6, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'D'  ,    600.0, 'd'     ,  400.0, 'B'      , 272.0);
                    PLC_A     = Bearing('name', 'PLC-A', 'type'  ,  'SRB', ...
                                        'a'   ,       a, 'C'     , 11164.0e3, 'e'      , 0.27, ...
                                        'X'   , [1.0 0.67], 'Y'  , [2.5 3.7], 'L_10'   , L_10, ...
                                        'K_x',    6.6e4, 'K_y'   ,  1.7e9, 'K_z'    , 1.1e9, ...
                                                         'K_beta',  5.6e5, 'K_gamma', 1.3e5, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'D'  ,   1030.0, 'd'     ,  710.0, 'B'      , 315.0);
                    PLC_B     = Bearing('name', 'PLC-B', 'type'  ,  'CRB', ...
                                        'a'   ,       a, 'C'     , 3690.0e3, 'e'      , 0.3, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x',    6.6e7, 'K_y'   ,  1.7e9, 'K_z'    , 1.1e9, ...
                                                         'K_beta',  5.6e5, 'K_gamma', 1.3e5, ...
                                        'C_x' ,      cx, 'C_y'   ,      cy, 'C_z'    ,     cz, ...
                                                         'C_beta',      cb, 'C_gamma',     cg, ...
                                        'D'  ,   1220.0, 'd'     , 1000.0, 'B'      , 128.0);
                    
                    brg = [PL_A,  PL_B,            ... % Planet
                           PLC_A, PLC_B,           ... % Carrier
                           Bearing('name', 'sun'), ... % Sun
                           Bearing('name', 'ring')];   % Ring
                           
                case 2
                    IMS_PL_A  = Bearing('name',  'IMS-PL-A', 'type'  , 'CRB' , ...
                                        'a'   ,       a, 'C'     , 2380.0e3, 'e'      , 0.15, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'D'   ,       520.0, 'd'     ,  380.0, 'B'      , 140.0);
                    IMS_PL_B  = Bearing('name',  'IMS-PL-B', 'type'  , 'CRB' , ...
                                        'a'   ,       a, 'C'     , 2380.0e3, 'e'      , 0.15, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'D'   ,       520.0, 'd'     ,  380.0, 'B'      , 140.0);
                    IMS_PLC_A = Bearing('name', 'IMS-PLC-A', 'type'  , 'CARB', ...
                                        'a'   ,       a, 'C'     , 8800.0e3, 'e'      , +inf, ...
                                        'X'   ,   [1.0 0.0], 'Y'  , [0.0 0.0], 'L_10'   , L_10, ...
                                        'K_x' ,       9.1e4, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'D'   ,      1030.0, 'd'     ,  710.0, 'B'      , 236.0);
                    IMS_PLC_B = Bearing('name', 'IMS-PLC-B', 'type'  , 'CRB' , ...
                                        'a'   ,       a, 'C'     , 4180.0e3, 'e'      , 0.3, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x' ,       9.1e7, 'K_y'   ,  6.0e7, 'K_z'    , 1.2e9, ...
                                                             'K_beta',  7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,          cx, 'C_y'   ,     cy, 'C_z'    ,    cz, ...
                                                             'C_beta',     cb, 'C_gamma',    cg, ...
                                        'D'   ,       870.0, 'd'     ,  600.0, 'B'      , 155.0);
                    
                    brg = [IMS_PL_A,  IMS_PL_B,    ... % Planet
                           IMS_PLC_A, IMS_PLC_B,   ... % Carrier
                           Bearing('name', 'sun'), ... % Sun
                           Bearing('name', 'ring')];   % Ring
                    
                case 3
                    IMS_A     = Bearing('name', 'IMS-A', 'type'  , 'CRB', ...
                                        'a'   ,       a, 'C'     , 1420.0e3, 'e'      , 0.15, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x' ,     0.0, 'K_y'   , 6.0e7, 'K_z'    , 1.2e9, ...
                                                         'K_beta', 7.5e4, 'K_gamma', 7.5e4, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'D'   ,   360.0, 'd'     , 200.0, 'B'      ,  98.0);
                    IMS_B     = Bearing('name', 'IMS-B', 'type'  , 'TRB', ...
                                        'a'   ,       a, 'C'     , 1467.0e3, 'e'      , 0.4, ...
                                        'X'   , [1.0 0.4], 'Y'  , [0.0 1.5], 'L_10'   , L_10, ...
                                        'K_x' ,   7.4e7, 'K_y'   , 5.0e8, 'K_z'    , 5.0e8, ...
                                                         'K_beta', 1.6e6, 'K_gamma', 1.8e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'D'   ,   460.0, 'd'     , 200.0, 'B'      , 100.0);
                    IMS_C     = Bearing('name', 'IMS-C', 'type'  , 'TRB', ...
                                        'a'   ,       a, 'C'     , 1467.0e3, 'e'      , 0.4, ...
                                        'X'   , [1.0 0.4], 'Y'  , [0.0 1.5], 'L_10'   , L_10, ...
                                        'K_x' ,   7.8e7, 'K_y'   , 7.4e8, 'K_z'    , 3.3e8, ...
                                                         'K_beta', 1.1e6, 'K_gamma', 2.5e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'D'   ,   460.0, 'd'     , 200.0, 'B'      , 100.0);
                    HS_A      = Bearing('name',  'HS-A', 'type'  , 'CRB', ...
                                        'a'   ,       a, 'C'     , 1460.0e3, 'e'      , 0.3, ...
                                        'X'   , [1.0 0.92], 'Y'  , [0.0 0.4], 'L_10'   , L_10, ...
                                        'K_x' ,   1.3e8, 'K_y'   , 8.2e8, 'K_z'    , 8.2e8, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'D'   ,   500.0, 'd'     , 400.0, 'B'      , 100.0);
                    HS_B      = Bearing('name',  'HS-B', 'type'  , 'TRB', ...
                                        'a'   ,       a, 'C'     , 1644.0e3, 'e'      , 0.43, ...
                                        'X'   , [1.0 0.4], 'Y'  , [0.0 1.4], 'L_10'   , L_10, ...
                                        'K_x' ,   6.7e7, 'K_y'   , 8.0e8, 'K_z'    , 1.3e8, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'D'   ,   550.0, 'd'     , 410.0, 'B'      ,  86.0);
                    HS_C      = Bearing('name',  'HS-C', 'type'  , 'TRB', ...
                                        'a'   ,       a, 'C'     , 1644.0e3, 'e'      , 0.43, ...
                                        'X'   , [1.0 0.4], 'Y'  , [0.0 1.4], 'L_10'   , L_10, ...
                                        'K_x' ,   8.0e7, 'K_y'   , 1.0e9, 'K_z'    , 7.3e7, ...
                                                         'K_beta', 1.7e5, 'K_gamma', 1.0e6, ...
                                        'C_x' ,      cx, 'C_y'   ,    cy, 'C_z'    ,    cz, ...
                                                         'C_beta',    cb, 'C_gamma',    cg, ...
                                        'D'   ,   550.0, 'd'     , 410.0, 'B'      ,  86.0);
                    
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
            
            sft = Shaft('d'      , d, ...
                        'L'      , L, ...
                        'bearing', bearing);
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
            u    =     97.0;
            m_r  =    110.0e3; % [kg]     , according to [2], Table 1-1
            J_g  =    534.116; % [kg-m^2] , according to [2], Table 5-1
            k_eq = 867637.0e3; % [N-m/rad], taken from [2], Table 5-1
            
            freq_free = 2.18; % [Hz], taken from [3], p. 45
            % Mean of the results obtained using FAST and ADAMS, respectively.
            freq_fix  = mean([0.6205 0.6094]); % [Hz], taken from [2], Table 9-1. 
            
            rho = Material().rho*1.0e9;
            
            J_r_fix   = k_eq/(2*pi*freq_fix)^2;
            J_r_free  = (k_eq*J_g*u^2)/(J_g*(2.0*pi*u*freq_free)^2 - k_eq);
            J_r_ratio = ((freq_free/freq_fix)^2 - 1.0)*J_g*u^2;
            
            R_cyl = sqrt(2.0*J_r_ratio/m_r); % [m], cylinder radius
            h_cyl = m_r/(rho*pi*R_cyl^2); % [m], cylinder height
            
            f = @(x)([(pi*rho/m_r)            .*x(1, :).*(1.0 - x(3, :).^2).*x(2, :)^2 - 1.0; ...
                      (pi*rho/(2.0*J_r_ratio)).*x(1, :).*(1.0 - x(3, :).^4).*x(2, :)^4 - 1.0]);
%             x0 = [2.0*h_cyl;
%                   0.5*R_cyl;
%                   0.5];
            x0 = rand(3,1);
            x = fmincon(@(x)(norm(f(x)))^2, x0, [], [], [], [], zeros(3, 1), [Inf, Inf, 1.0]');
            r_ext = x(2);           r_int = x(3)*r_ext;           h = x(1);
            
            J_other = (1.0/12.0)*m_r*(3.0*(r_ext^2 + r_int^2) + h^2);

            fprintf('Rotor mass moment of inertia using:\n');
            fprintf('\t Rigid free-fixed resonance: %3.4e [Hz]\t %3.4e [kg-m^2]\n', freq_fix, J_r_fix);
            fprintf('\t Rigid free-free  resonance: %3.4e [Hz]\t %3.4e [kg-m^2]\n', freq_free, J_r_free);
            fprintf('\t Ratio of the frequencies above: \t\t %3.4e [kg-m^2]\n', J_r_ratio);
            fprintf('Rotor mass moment of inertia on the other axes: %3.4e [kg-m^2]\n', J_other);
            fprintf('Rotor dimensions assuming:\n');
            fprintf('\t Cylindric geometry: R = %.3f \t h = %.3f [m]\n', R_cyl, h_cyl);
            fprintf('\t Cylindric tube geometry, opt.: R_ext = %.3f \t R_int = %.3f \t h = %.3f [m]\n', r_ext, r_int, h);
        end
    end
    
    methods
%         function update_subvar(obj, varargin)
%             default = {'simulation_time'     , 200.0, ... % [s]
%                        'inertia_flag'        , true, ...
%                        'base_excitation_flag', true, ...
%                        'generator_flag'      , 2, ... % PI control
%                        'mesh_flag'           , 225};
%             
%             default = scaling_factor.process_varargin(default, varargin);
%             
%             gamma_module = mean([obj.gamma('m_n1') obj.gamma('m_n2') obj.gamma('m_n3')]);
%             gamma_load = power(gamma_module, 2.0);
%             gamma_Torque = power(gamma_module, 3.0);
%             
%             k_SP = zeros(2, 1);
%             k_RP = k_SP;
%             for idx = 1:2
%                 k_SP(idx) = obj.stage(idx).sub_set('sun_planet').k_mesh;
%                 k_RP(idx) = obj.stage(idx).sub_set('planet_ring').k_mesh;
%             end
%             
%             k_PW = obj.stage(idx + 1).k_mesh;
% 
%             new_ID = fopen('@NREL_5MW\NREL_5MW.scaled.subvar', 'w');
%             
%             fprintf(new_ID, "!file.version=3.5! Removing this line will make the file unreadable\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "!**********************************************************************\n");
%             fprintf(new_ID, "! SubVars\n");
%             fprintf(new_ID, "!**********************************************************************\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_loading                  )                           ! $SVG_loading\n");
%             fprintf(new_ID, "   subvar($_gamma_Power, str= '%.5e')                                                          ! $SVG_loading.$_gamma_Power\n", obj.gamma('P'));
%             fprintf(new_ID, "   subvar($_gamma_Speed, str= '1.0')                                                          ! $SVG_loading.$_gamma_Speed\n");
%             fprintf(new_ID, "   subvar($_gamma_Torque, str= '$SVG_loading.$_gamma_Power/$SVG_loading.$_gamma_Speed')       ! $SVG_loading.$_gamma_Torque\n");
%             fprintf(new_ID, "   subvar($_gamma_Force, str= 'pow($SVG_loading.$_gamma_Power, 2.0/3.0)/$SVG_loading.$_gamma_Speed') ! $SVG_loading.$_gamma_Force\n");
%             fprintf(new_ID, "   subvar($_rated_power_ref, str= '5.0e3 kW')                                                 ! $SVG_loading.$_rated_power_ref\n");
%             fprintf(new_ID, "   subvar($_rated_power, str= '$SVG_loading.$_rated_power_ref*$SVG_loading.$_gamma_Power')    ! $SVG_loading.$_rated_power\n");
%             fprintf(new_ID, "   subvar($_inertia_flag, str= '%d', discr.desc (   1) = 'yes (value)', discr.desc (   2) = 'no (unit)', discr.str (   1) = '1', discr.str (   2) = '0') ! $SVG_loading.$_inertia_flag\n", default.inertia_flag);
%             fprintf(new_ID, "   subvar($_base_excitation_flag, str= '%d', discr.desc (   1) = 'ON', discr.desc (   2) = 'OFF', discr.str (   1) = '1', discr.str (   2) = '0') ! $SVG_loading.$_base_excitation_flag\n", default.base_excitation_flag);
%             fprintf(new_ID, "   subvar($_generator_flag, str= '%d', discr.desc (   1) = 'Prop. control', discr.desc (   2) = 'PI control', discr.desc (   3) = 'Time series, omega_gen', discr.str (   1) = '1', discr.str (   2) = '2', discr.str (   3) = '3') ! $SVG_loading.$_generator_flag\n", default.generator_flag);
%             fprintf(new_ID, "   subvar($_mesh_flag, str= '%d', discr.desc (   1) = 'advanced', discr.desc (   2) = 'basic', discr.str (   1) = '225', discr.str (   2) = '204') ! $SVG_loading.$_mesh_flag\n", default.mesh_flag);
%             fprintf(new_ID, "   subvar($_simulation_time, str= '%.5e s')                                                   ! $SVG_loading.$_simulation_time\n", default.simulation_time);
%             fprintf(new_ID, "   subvar($_gamma_Disp, str= '%.5e')                                                           ! $SVG_loading.$_gamma_Disp\n", obj.gamma('m_n1'));
%             fprintf(new_ID, "subvargroup.end (                   $SVG_loading                  )                           ! $SVG_loading\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_rotor                    )                           ! $SVG_rotor\n");
%             fprintf(new_ID, "   subvar($_gamma_M, str= '1.0')                                                              ! $SVG_rotor.$_gamma_M\n");
%             fprintf(new_ID, "   subvar($_gamma_J, str= '%.5e')                                                                ! $SVG_rotor.$_gamma_J\n", obj.gamma('J_R'));
%             fprintf(new_ID, "   subvar($_mass_ref, str= 'IF($SVG_loading.$_inertia_flag == 1)\n{\n110.0e3 kg\n}\nELSE\n{\n1.0 kg\n}') ! $SVG_rotor.$_mass_ref\n");
%             fprintf(new_ID, "   subvar($_mass_mom_inertia_ref, str= 'IF($SVG_loading.$_inertia_flag== 1)\n{\n57231535.0 kg m^2\n}\nELSE\n{\n1.0 kg m^2\n}') ! $SVG_rotor.$_mass_mom_inertia_ref\n");
%             fprintf(new_ID, "   subvar($_length_ref, str= '250.0 mm')                                                      ! $SVG_rotor.$_length_ref\n");
%             fprintf(new_ID, "   subvar($_diameter_ref, str= '840.0 mm')                                                    ! $SVG_rotor.$_diameter_ref\n");
%             fprintf(new_ID, "   subvar($_mass, str= '$SVG_rotor.$_mass_ref*$SVG_rotor.$_gamma_M')                          ! $SVG_rotor.$_mass\n");
%             fprintf(new_ID, "   subvar($_mass_mom_inertia, str= '$SVG_rotor.$_mass_mom_inertia_ref*$SVG_rotor.$_gamma_J')  ! $SVG_rotor.$_mass_mom_inertia\n");
%             fprintf(new_ID, "   subvar($_length, str= '$SVG_rotor.$_length_ref*$SVG_main_shaft.$_gamma_L')                 ! $SVG_rotor.$_length\n");
%             fprintf(new_ID, "   subvar($_diameter, str= '$SVG_rotor.$_diameter_ref*$SVG_main_shaft.$_gamma_d')             ! $SVG_rotor.$_diameter\n");
%             fprintf(new_ID, "   subvar($_rated_speed, str= '12.1 rpm')                                                     ! $SVG_rotor.$_rated_speed\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_rotor                    )                           ! $SVG_rotor\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_main_shaft               )                           ! $SVG_main_shaft\n");
%             fprintf(new_ID, "   subvar($_gamma_L, str= '%.5e')                                                                ! $SVG_main_shaft.$_gamma_L\n", obj.gamma('L_S'));
%             fprintf(new_ID, "   subvar($_gamma_d, str= '%.5e')                                                                ! $SVG_main_shaft.$_gamma_d\n", obj.gamma('d_S'));
%             fprintf(new_ID, "   subvar($_length_ref, str= '2000.0 mm')                                                     ! $SVG_main_shaft.$_length_ref\n");
%             fprintf(new_ID, "   subvar($_diameter_ref, str= '700.0 mm')                                                    ! $SVG_main_shaft.$_diameter_ref\n");
%             fprintf(new_ID, "   subvar($_x_A_ref, str= '333.0 mm')                                                         ! $SVG_main_shaft.$_x_A_ref\n");
%             fprintf(new_ID, "   subvar($_x_B_ref, str= '1166.0 mm')                                                        ! $SVG_main_shaft.$_x_B_ref\n");
%             fprintf(new_ID, "   subvar($_length, str= '$SVG_main_shaft.$_length_ref*$SVG_main_shaft.$_gamma_L')            ! $SVG_main_shaft.$_length\n");
%             fprintf(new_ID, "   subvar($_diameter, str= '$SVG_main_shaft.$_diameter_ref*$SVG_main_shaft.$_gamma_d')        ! $SVG_main_shaft.$_diameter\n");
%             fprintf(new_ID, "   subvar($_x_A, str= '$SVG_main_shaft.$_x_A_ref*$SVG_main_shaft.$_gamma_L')                  ! $SVG_main_shaft.$_x_A\n");
%             fprintf(new_ID, "   subvar($_x_B, str= '$SVG_main_shaft.$_x_B_ref*$SVG_main_shaft.$_gamma_L')                  ! $SVG_main_shaft.$_x_B\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_main_shaft               )                           ! $SVG_main_shaft\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_stage_01                 )                           ! $SVG_stage_01\n");
%             fprintf(new_ID, "   subvar($_gamma_mn, str= '%.5e')                                                               ! $SVG_stage_01.$_gamma_mn\n", obj.gamma('m_n1'));
%             fprintf(new_ID, "   subvar($_gamma_b, str= '%.5e')                                                                ! $SVG_stage_01.$_gamma_b\n", obj.gamma('b_1'));
%             fprintf(new_ID, "   subvar($_gamma_L, str= '%.5e')                                                                ! $SVG_stage_01.$_gamma_L\n", obj.gamma('L_1'));
%             fprintf(new_ID, "   subvar($_gamma_d, str= '%.5e')                                                                ! $SVG_stage_01.$_gamma_d\n", obj.gamma('d_1'));
%             fprintf(new_ID, "   subvar($_k_SP, str= '%.5e N/m')                                                           ! $SVG_stage_01.$_k_SP\n", k_SP(1));
%             fprintf(new_ID, "   subvar($_k_RP, str= '%.5e N/m')                                                           ! $SVG_stage_01.$_k_RP\n", k_RP(1));
%             fprintf(new_ID, "   subvar($_normal_module_ref, str= '45 mm')                                                  ! $SVG_stage_01.$_normal_module_ref\n");
%             fprintf(new_ID, "   subvar($_center_distance_ref, str= '863.0 mm')                                             ! $SVG_stage_01.$_center_distance_ref\n");
%             fprintf(new_ID, "   subvar($_normal_backlash_ref, str= '0.40 mm')                                              ! $SVG_stage_01.$_normal_backlash_ref\n");
% %             fprintf(new_ID, "   subvar($_normal_module, str= '$SVG_stage_01.$_normal_module_ref*$SVG_stage_01.$_gamma_mn') ! $SVG_stage_01.$_normal_module\n");
% %             fprintf(new_ID, "   subvar($_center_distance, str= '$SVG_stage_01.$_center_distance_ref*$SVG_stage_01.$_gamma_mn') ! $SVG_stage_01.$_center_distance\n");
%             fprintf(new_ID, "   subvar($_normal_module, str= '%.5e mm') ! $SVG_stage_01.$_normal_module\n", obj.stage(1).m_n);
%             fprintf(new_ID, "   subvar($_center_distance, str= '%.5e mm') ! $SVG_stage_01.$_center_distance\n", obj.stage(1).a_w);
%             fprintf(new_ID, "   subvar($_normal_backlash, str= '$SVG_stage_01.$_normal_backlash_ref*$SVG_stage_01.$_gamma_mn') ! $SVG_stage_01.$_normal_backlash\n");
%             fprintf(new_ID, "   subvar($_length, str= '1.8*$SVG_stage_01.$SVG_planet.$_flank_width + $SVG_stage_01.$SVG_shaft.$_length') ! $SVG_stage_01.$_length\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_sun                      )                           ! $SVG_stage_01.$SVG_sun\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '491.0 mm')                                              ! $SVG_stage_01.$SVG_sun.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_bore_diameter_ref, str= '400.0 mm')                                            ! $SVG_stage_01.$SVG_sun.$_bore_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_01.$SVG_sun.$_flank_width_ref*$SVG_stage_01.$_gamma_b') ! $SVG_stage_01.$SVG_sun.$_flank_width\n");
%             fprintf(new_ID, "      subvar($_bore_diameter, str= '$SVG_stage_01.$SVG_sun.$_bore_diameter_ref*$SVG_stage_01.$_gamma_mn') ! $SVG_stage_01.$SVG_sun.$_bore_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_sun                      )                           ! $SVG_stage_01.$SVG_sun\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_planet                   )                           ! $SVG_stage_01.$SVG_planet\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '$SVG_stage_01.$SVG_sun.$_flank_width_ref')              ! $SVG_stage_01.$SVG_planet.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_bore_diameter_ref, str= '400.0 mm')                                            ! $SVG_stage_01.$SVG_planet.$_bore_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_01.$SVG_planet.$_flank_width_ref*$SVG_stage_01.$_gamma_b') ! $SVG_stage_01.$SVG_planet.$_flank_width\n");
%             fprintf(new_ID, "      subvar($_bore_diameter, str= '$SVG_stage_01.$SVG_planet.$_bore_diameter_ref*$SVG_stage_01.$_gamma_mn') ! $SVG_stage_01.$SVG_planet.$_bore_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_planet                   )                           ! $SVG_stage_01.$SVG_planet\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_ring                     )                           ! $SVG_stage_01.$SVG_ring\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '$SVG_stage_01.$SVG_sun.$_flank_width_ref')              ! $SVG_stage_01.$SVG_ring.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_01.$SVG_ring.$_flank_width_ref*$SVG_stage_01.$_gamma_b') ! $SVG_stage_01.$SVG_ring.$_flank_width\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_ring                     )                           ! $SVG_stage_01.$SVG_ring\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_shaft                    )                           ! $SVG_stage_01.$SVG_shaft\n");
%             fprintf(new_ID, "      subvar($_length_ref, str= '500.0 mm')                                                   ! $SVG_stage_01.$SVG_shaft.$_length_ref\n");
%             fprintf(new_ID, "      subvar($_diameter_ref, str= '533.0 mm')                                                 ! $SVG_stage_01.$SVG_shaft.$_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_length, str= '$SVG_stage_01.$SVG_shaft.$_length_ref*$SVG_stage_01.$_gamma_L')  ! $SVG_stage_01.$SVG_shaft.$_length\n");
%             fprintf(new_ID, "      subvar($_diameter, str= '$SVG_stage_01.$SVG_shaft.$_diameter_ref*$SVG_stage_01.$_gamma_d') ! $SVG_stage_01.$SVG_shaft.$_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_shaft                    )                           ! $SVG_stage_01.$SVG_shaft\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_stage_01                 )                           ! $SVG_stage_01\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_stage_02                 )                           ! $SVG_stage_02\n");
%             fprintf(new_ID, "   subvar($_gamma_mn, str= '%.5e')                                                               ! $SVG_stage_02.$_gamma_mn\n", obj.gamma('m_n2'));
%             fprintf(new_ID, "   subvar($_gamma_b, str= '%.5e')                                                                ! $SVG_stage_02.$_gamma_b\n", obj.gamma('b_2'));
%             fprintf(new_ID, "   subvar($_gamma_L, str= '%.5e')                                                                ! $SVG_stage_02.$_gamma_L\n", obj.gamma('L_2'));
%             fprintf(new_ID, "   subvar($_gamma_d, str= '%.5e')                                                                ! $SVG_stage_02.$_gamma_d\n", obj.gamma('d_2'));
%             fprintf(new_ID, "   subvar($_k_SP, str= '%.5e N/m')                                                           ! $SVG_stage_02.$_k_SP\n", k_SP(2));
%             fprintf(new_ID, "   subvar($_k_RP, str= '%.5e N/m')                                                           ! $SVG_stage_02.$_k_RP\n", k_RP(2));
%             fprintf(new_ID, "   subvar($_normal_module_ref, str= '21.0 mm')                                                ! $SVG_stage_02.$_normal_module_ref\n");
%             fprintf(new_ID, "   subvar($_center_distance_ref, str= '584.0 mm')                                             ! $SVG_stage_02.$_center_distance_ref\n");
%             fprintf(new_ID, "   subvar($_normal_backlash_ref, str= '0.40 mm')                                              ! $SVG_stage_02.$_normal_backlash_ref\n");
% %             fprintf(new_ID, "   subvar($_normal_module, str= '$SVG_stage_02.$_normal_module_ref*$SVG_stage_02.$_gamma_mn') ! $SVG_stage_02.$_normal_module\n");
% %             fprintf(new_ID, "   subvar($_center_distance, str= '$SVG_stage_02.$_center_distance_ref*$SVG_stage_02.$_gamma_mn') ! $SVG_stage_02.$_center_distance\n");
%             fprintf(new_ID, "   subvar($_normal_module, str= '%.5e mm') ! $SVG_stage_02.$_normal_module\n", obj.stage(2).m_n);
%             fprintf(new_ID, "   subvar($_center_distance, str= '%.5e mm') ! $SVG_stage_02.$_center_distance\n", obj.stage(2).a_w);
%             fprintf(new_ID, "   subvar($_normal_backlash, str= '$SVG_stage_02.$_normal_backlash_ref*$SVG_stage_02.$_gamma_mn') ! $SVG_stage_02.$_normal_backlash\n");
%             fprintf(new_ID, "   subvar($_length, str= '1.8*$SVG_stage_02.$SVG_planet.$_flank_width + $SVG_stage_02.$SVG_shaft.$_length') ! $SVG_stage_02.$_length\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_sun                      )                           ! $SVG_stage_02.$SVG_sun\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '550.0 mm')                                              ! $SVG_stage_02.$SVG_sun.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_bore_diameter_ref, str= '200.0 mm')                                            ! $SVG_stage_02.$SVG_sun.$_bore_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_02.$SVG_sun.$_flank_width_ref*$SVG_stage_02.$_gamma_b') ! $SVG_stage_02.$SVG_sun.$_flank_width\n");
%             fprintf(new_ID, "      subvar($_bore_diameter, str= '$SVG_stage_02.$SVG_sun.$_bore_diameter_ref*$SVG_stage_02.$_gamma_mn') ! $SVG_stage_02.$SVG_sun.$_bore_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_sun                      )                           ! $SVG_stage_02.$SVG_sun\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_planet                   )                           ! $SVG_stage_02.$SVG_planet\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '$SVG_stage_02.$SVG_sun.$_flank_width_ref')              ! $SVG_stage_02.$SVG_planet.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_bore_diameter_ref, str= '380.0 mm')                                            ! $SVG_stage_02.$SVG_planet.$_bore_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_02.$SVG_planet.$_flank_width_ref*$SVG_stage_02.$_gamma_b') ! $SVG_stage_02.$SVG_planet.$_flank_width\n");
%             fprintf(new_ID, "      subvar($_bore_diameter, str= '$SVG_stage_02.$SVG_planet.$_bore_diameter_ref*$SVG_stage_02.$_gamma_mn') ! $SVG_stage_02.$SVG_planet.$_bore_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_planet                   )                           ! $SVG_stage_02.$SVG_planet\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_ring                     )                           ! $SVG_stage_02.$SVG_ring\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '$SVG_stage_02.$SVG_sun.$_flank_width_ref')              ! $SVG_stage_02.$SVG_ring.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_02.$SVG_ring.$_flank_width_ref*$SVG_stage_02.$_gamma_b') ! $SVG_stage_02.$SVG_ring.$_flank_width\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_ring                     )                           ! $SVG_stage_02.$SVG_ring\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_shaft                    )                           ! $SVG_stage_02.$SVG_shaft\n");
%             fprintf(new_ID, "      subvar($_length_ref, str= '666.0 mm')                                                   ! $SVG_stage_02.$SVG_shaft.$_length_ref\n");
%             fprintf(new_ID, "      subvar($_diameter_ref, str= '333.0 mm')                                                 ! $SVG_stage_02.$SVG_shaft.$_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_length, str= '$SVG_stage_02.$SVG_shaft.$_length_ref*$SVG_stage_02.$_gamma_L')  ! $SVG_stage_02.$SVG_shaft.$_length\n");
%             fprintf(new_ID, "      subvar($_diameter, str= '$SVG_stage_02.$SVG_shaft.$_diameter_ref*$SVG_stage_02.$_gamma_d') ! $SVG_stage_02.$SVG_shaft.$_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_shaft                    )                           ! $SVG_stage_02.$SVG_shaft\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_stage_02                 )                           ! $SVG_stage_02\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_stage_03                 )                           ! $SVG_stage_03\n");
%             fprintf(new_ID, "   subvar($_gamma_mn, str= '%.5e')                                                               ! $SVG_stage_03.$_gamma_mn\n", obj.gamma('m_n3'));
%             fprintf(new_ID, "   subvar($_gamma_b, str= '%.5e')                                                                ! $SVG_stage_03.$_gamma_b\n", obj.gamma('b_3'));
%             fprintf(new_ID, "   subvar($_gamma_L, str= '%.5e')                                                                ! $SVG_stage_03.$_gamma_L\n", obj.gamma('L_3'));
%             fprintf(new_ID, "   subvar($_gamma_d, str= '%.5e')                                                                ! $SVG_stage_03.$_gamma_d\n", obj.gamma('d_3'));
%             fprintf(new_ID, "   subvar($_k_PW, str= '%.5e N/m')                                                           ! $SVG_stage_03.$_k_PW\n", k_PW);
%             fprintf(new_ID, "   subvar($_normal_module_ref, str= '14 mm')                                                  ! $SVG_stage_03.$_normal_module_ref\n");
%             fprintf(new_ID, "   subvar($_center_distance_ref, str= '861.0 mm')                                             ! $SVG_stage_03.$_center_distance_ref\n");
%             fprintf(new_ID, "   subvar($_normal_backlash_ref, str= '0.40 mm')                                              ! $SVG_stage_03.$_normal_backlash_ref\n");
% %             fprintf(new_ID, "   subvar($_normal_module, str= '$SVG_stage_03.$_normal_module_ref*$SVG_stage_03.$_gamma_mn') ! $SVG_stage_03.$_normal_module\n");
% %             fprintf(new_ID, "   subvar($_center_distance, str= '$SVG_stage_03.$_center_distance_ref*$SVG_stage_03.$_gamma_mn') ! $SVG_stage_03.$_center_distance\n");
%             fprintf(new_ID, "   subvar($_normal_module, str= '%.5e mm') ! $SVG_stage_03.$_normal_module\n", obj.stage(3).m_n);
%             fprintf(new_ID, "   subvar($_center_distance, str= '%.5e mm') ! $SVG_stage_03.$_center_distance\n", obj.stage(3).a_w);
%             fprintf(new_ID, "   subvar($_normal_backlash, str= '$SVG_stage_03.$_normal_backlash_ref*$SVG_stage_03.$_gamma_mn') ! $SVG_stage_03.$_normal_backlash\n");
%             fprintf(new_ID, "   subvar($_length, str= '1.8*$SVG_stage_03.$SVG_pinion.$_flank_width + $SVG_stage_03.$SVG_shaft.$_length') ! $SVG_stage_03.$_length\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_pinion                   )                           ! $SVG_stage_03.$SVG_pinion\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '360.0 mm')                                              ! $SVG_stage_03.$SVG_pinion.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_bore_diameter_ref, str= '200.0 mm')                                            ! $SVG_stage_03.$SVG_pinion.$_bore_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_03.$SVG_pinion.$_flank_width_ref*$SVG_stage_03.$_gamma_b') ! $SVG_stage_03.$SVG_pinion.$_flank_width\n");
%             fprintf(new_ID, "      subvar($_bore_diameter, str= '$SVG_stage_03.$SVG_pinion.$_bore_diameter_ref*$SVG_stage_03.$_gamma_mn') ! $SVG_stage_03.$SVG_pinion.$_bore_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_pinion                   )                           ! $SVG_stage_03.$SVG_pinion\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_wheel                    )                           ! $SVG_stage_03.$SVG_wheel\n");
%             fprintf(new_ID, "      subvar($_flank_width_ref, str= '$SVG_stage_03.$SVG_pinion.$_flank_width_ref')           ! $SVG_stage_03.$SVG_wheel.$_flank_width_ref\n");
%             fprintf(new_ID, "      subvar($_bore_diameter_ref, str= '400.0 mm')                                            ! $SVG_stage_03.$SVG_wheel.$_bore_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_flank_width, str= '$SVG_stage_03.$SVG_wheel.$_flank_width_ref*$SVG_stage_03.$_gamma_b') ! $SVG_stage_03.$SVG_wheel.$_flank_width\n");
%             fprintf(new_ID, "      subvar($_bore_diameter, str= '$SVG_stage_03.$SVG_wheel.$_bore_diameter_ref*$SVG_stage_03.$_gamma_mn') ! $SVG_stage_03.$SVG_wheel.$_bore_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_wheel                    )                           ! $SVG_stage_03.$SVG_wheel\n");
%             fprintf(new_ID, "   subvargroup.begin (              $SVG_shaft                    )                           ! $SVG_stage_03.$SVG_shaft\n");
%             fprintf(new_ID, "      subvar($_length_ref, str= '1000.0 mm')                                                  ! $SVG_stage_03.$SVG_shaft.$_length_ref\n");
%             fprintf(new_ID, "      subvar($_diameter_ref, str= '333.0 mm')                                                 ! $SVG_stage_03.$SVG_shaft.$_diameter_ref\n");
%             fprintf(new_ID, "      subvar($_length, str= '$SVG_stage_03.$SVG_shaft.$_length_ref*$SVG_stage_03.$_gamma_L')  ! $SVG_stage_03.$SVG_shaft.$_length\n");
%             fprintf(new_ID, "      subvar($_diameter, str= '$SVG_stage_03.$SVG_shaft.$_diameter_ref*$SVG_stage_03.$_gamma_d') ! $SVG_stage_03.$SVG_shaft.$_diameter\n");
%             fprintf(new_ID, "   subvargroup.end (                $SVG_shaft                    )                           ! $SVG_stage_03.$SVG_shaft\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_stage_03                 )                           ! $SVG_stage_03\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_generator                )                           ! $SVG_generator\n");
%             fprintf(new_ID, "   subvar($_gamma_M, str= '1.0')                                                              ! $SVG_generator.$_gamma_M\n");
%             fprintf(new_ID, "   subvar($_gamma_J, str= '1.0')                                                              ! $SVG_generator.$_gamma_J\n");
%             fprintf(new_ID, "   subvar($_mass_ref, str= 'IF($SVG_loading.$_inertia_flag == 1)\n{\n1900.0 kg\n}\nELSE\n{\n1.0 kg\n}') ! $SVG_generator.$_mass_ref\n");
%             fprintf(new_ID, "   subvar($_mass_mom_inertia_ref, str= 'IF($SVG_loading.$_inertia_flag  == 1)\n{\n534.116 kg m^2\n}\nELSE\n{\n1.0 kg m^2\n}') ! $SVG_generator.$_mass_mom_inertia_ref\n");
%             fprintf(new_ID, "   subvar($_length_ref, str= '$SVG_rotor.$_length_ref')                                       ! $SVG_generator.$_length_ref\n");
%             fprintf(new_ID, "   subvar($_diameter_ref, str= '$SVG_rotor.$_diameter_ref')                                   ! $SVG_generator.$_diameter_ref\n");
%             fprintf(new_ID, "   subvar($_K_P_ref, str= '2200.0', desc (   1) = 'Proportional gain')                        ! $SVG_generator.$_K_P_ref\n");
%             fprintf(new_ID, "   subvar($_K_I_ref, str= '220.0', desc (   1) = 'Integral gain')                             ! $SVG_generator.$_K_I_ref\n");
%             fprintf(new_ID, "   subvar($_mass, str= '$SVG_generator.$_mass_ref*$SVG_generator.$_gamma_M')                  ! $SVG_generator.$_mass\n");
%             fprintf(new_ID, "   subvar($_mass_mom_inertia, str= '$SVG_generator.$_mass_mom_inertia_ref*$SVG_generator.$_gamma_J') ! $SVG_generator.$_mass_mom_inertia\n");
%             fprintf(new_ID, "   subvar($_length, str= '$SVG_generator.$_length_ref*$SVG_stage_03.$_gamma_L')               ! $SVG_generator.$_length\n");
%             fprintf(new_ID, "   subvar($_diameter, str= '$SVG_generator.$_diameter_ref*$SVG_main_shaft.$_gamma_d')         ! $SVG_generator.$_diameter\n");
%             fprintf(new_ID, "   subvar($_K_P, str= '$SVG_generator.$_K_P_ref*$SVG_loading.$_gamma_Torque', desc (   1) = 'Proportional gain') ! $SVG_generator.$_K_P\n");
%             fprintf(new_ID, "   subvar($_K_I, str= '$SVG_generator.$_K_I_ref*$SVG_loading.$_gamma_Torque', desc (   1) = 'Integral gain') ! $SVG_generator.$_K_I\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_generator                )                           ! $SVG_generator\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_housing                  )                           ! $SVG_housing\n");
%             fprintf(new_ID, "   subvar($_length, str= '$SVG_stage_01.$_length + $SVG_stage_02.$_length + $SVG_stage_03.$_length') ! $SVG_housing.$_length\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_housing                  )                           ! $SVG_housing\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "subvargroup.begin (                 $SVG_bed_plate                )                           ! $SVG_bed_plate\n");
%             fprintf(new_ID, "   subvar($_length_ref, str= '8000.0 mm')                                                     ! $SVG_bed_plate.$_length_ref\n");
%             fprintf(new_ID, "   subvar($_width_ref, str= '4000.0 mm')                                                      ! $SVG_bed_plate.$_width_ref\n");
%             fprintf(new_ID, "   subvar($_length, str= '1.75*$SVG_housing.$_length')                                        ! $SVG_bed_plate.$_length\n");
%             fprintf(new_ID, "   subvar($_width, str= '$SVG_bed_plate.$_width_ref*$SVG_stage_01.$_gamma_mn')                ! $SVG_bed_plate.$_width\n");
%             fprintf(new_ID, "subvargroup.end (                   $SVG_bed_plate                )                           ! $SVG_bed_plate\n");
%             fprintf(new_ID, "\n");
%             fprintf(new_ID, "                                                                                                                                                                                                                                                                                                                                                                                                    \n");
%             
%             fclose(new_ID);
%             
%         end
%         
    end
    
end
