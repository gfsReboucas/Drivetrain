classdef WTDB_5MW < Drivetrain
    
    properties
        gamma_P (1, 1) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0; % Scaling factor for the rated power
        gamma_n (1, 1) {mustBeNumeric, mustBeFinite, mustBePositive} = 1.0; % Scaling factor for the rotor speed
    end
    
    methods
        function obj = WTDB_5MW(varargin)
            gamma = {'P'   , 1.0, 'n'   , 1.0, ...
                     'm_n1', 1.0, 'b_1' , 1.0, 'd_1' , 1.0, 'L_1' , 1.0, ...
                     'm_n2', 1.0, 'b_2' , 1.0, 'd_2' , 1.0, 'L_2' , 1.0, ...
                     'm_n3', 1.0, 'b_3' , 1.0, 'd_3' , 1.0, 'L_3' , 1.0, ...
                     'J_R' , 1.0, 'J_G' , 1.0, ...
                     'd_S' , 1.0, 'L_S' , 1.0};
            
            gamma = process_varargin(gamma, varargin);
            gamma = scaling_factor(gamma);
                        
            N_st = 3;
            stage = repmat(Gear_Set(), 1, N_st);
            
            for idx = 1:N_st
                stg = WTDB_5MW.gear_stage(idx);
                
                gamma_idx = gamma.ends_with(num2str(idx));
                stage(idx) = stg.scale_by(gamma_idx);
            end
            
            P_r = gamma('P')* 5.0e3; % [kW], Rated power
            n_r = gamma('n')*12.1;   % [1/min.], Input speed
            
            LSS = WTDB_5MW.shaft(0);
            
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
                           'dynamic_model', @Kahraman_94);
%                            'dynamic_model', @Lin_Parker_99);

            obj.gamma = scaling_factor(gamma.name, gamma.value);
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
                    p    =    4;         % [-],    Number of planet gears
                    m_n  =   25.0;       % [mm],   Normal module
                    beta =    7.0;       % [deg.], Helix angle (at reference cylinder)
                    b    = [540.0, ...
                            520.0, ...
                            520.0];     % [mm],   Face width
                    a_w  =  724.995;       % [mm],   Center distance
                    z_s  =   23;         % [-],    Number of teeth (sun)    [WHEEL]
                    z_p  =   34;         % [-],    Number of teeth (planet) [PINION]
                    z_r  =   93;         % [-],    Number of teeth (ring)
                    x_s  =    0.14;      % [-],    Profile shift coefficient (sun)
                    x_p  =    0.156;     % [-],    Profile shift coefficient (planet)
                    x_r  =    0.4937;    % [-],    Profile shift coefficient (ring)
                    k_s  =  -10.861/m_n; % [-],    Tip alteration coefficient (sun)
                    k_p  =  -10.861/m_n; % [-],    Tip alteration coefficient (planet)
                    k_r  =    0.0;       % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs = 200.0/(z_s*m_n./cosd(beta));
                    bore_Rp = 430.0/(z_p*m_n./cosd(beta));
                    bore_Rr =  1.2;
                    
                    z = [z_s z_p z_r];
                    x = [x_s x_p x_r];
                    k = [k_s k_p k_r];
                    bore_ratio = [bore_Rs bore_Rp bore_Rr];
                    config = 'planetary';
                    
                case 2
                    p    =   3;        % [-],    Number of planet gears
                    m_n  =  18.0;      % [mm],   Normal module
                    beta =   7.0;      % [deg.], Helix angle (at reference cylinder)
                    b    = [410.0, ...
                            390.0, ...
                            390.0];    % [mm],   Face width
                    a_w  = 584.0;      % [mm],   Center distance
                    z_s  =  20;        % [-],    Number of teeth (sun)    [PINION]
                    z_p  =  37;        % [-],    Number of teeth (planet) [WHEEL]
                    z_r  =  93;        % [-],    Number of teeth (ring)
                    x_s  =   0.42;     % [-],    Profile shift coefficient (sun)
                    x_p  =   0.3741;   % [-],    Profile shift coefficient (planet)
                    x_r  =   0.117;    % [-],    Profile shift coefficient (ring)
                    k_s  =  -1.75/m_n; % [-],    Tip alteration coefficient (sun)
                    k_p  =  -1.75/m_n; % [-],    Tip alteration coefficient (planet)
                    k_r  =   0.0;      % [-],    Tip alteration coefficient (ring)
                    
                    bore_Rs = 200.0/(z_s*m_n./cosd(beta));
                    bore_Rp = 430.0/(z_p*m_n./cosd(beta));
                    bore_Rr =   1.2;

                    z = [z_s z_p z_r];
                    x = [x_s x_p x_r];
                    k = [k_s k_p k_r];
                    bore_ratio = [bore_Rs bore_Rp bore_Rr];
                    config = 'planetary';

                case 3
                    p    =   1;         % [-],    Number of planet gears
                    m_n  =  16.0;       % [mm],   Normal module
                    beta =   7.0;       % [deg.], Helix angle (at reference cylinder)
                    b    = [300.0, ...
                            280.0];     % [mm],   Face width
                    a_w  = 861.0;       % [mm],   Center distance
                    z_1  =  19;         % [-],    Number of teeth (pinion)
                    z_2  =  63;         % [-],    Number of teeth (wheel)
                    x_1  =   0.52;      % [-],    Profile shift coefficient (pinion)
                    x_2  =   0.789;     % [-],    Profile shift coefficient (wheel)
                    k_1  =  -0.938/m_n; % [-],    Tip alteration coefficient (pinion)
                    k_2  =  -0.938/m_n; % [-],    Tip alteration coefficient (wheel)
                    
                    bore_R1 = 100.0/(z_1*m_n./cosd(beta));
                    bore_R2 = 520.0/(z_2*m_n./cosd(beta));

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
        
    end
    
end
