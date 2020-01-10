classdef Bearing
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
    
    properties
%         X;
%         Y;
%         a;
        name(1, :) string;    % [-],       Bearing designation
        type(1, :) string {mustBeMember(type, ["CARB", "CRB", "SRB", "TRB", "none"])} = "none";
                                              % CARB: ?
                                                      % CRB: Cylindrical Roller Bearing
                                                             % SRB: Spherical Roller Bearing
                                                                    % TRB: Taper Roller Bearing

        K_x;     % [N/m],     Translational stiffness, x axis
        K_y;     % [N/m],     Translational stiffness, y axis
        K_z;     % [N/m],     Translational stiffness, z axis
        K_alpha; % [N-m/rad], Torsional stiffness, x axis (rot.)
        K_beta;  % [N-m/rad], Torsional stiffness, y axis
        K_gamma; % [N-m/rad], Torsional stiffness, z axis

        OD;      % [mm],      Outer diameter
        ID;      % [mm],      Inner diameter
        B;       % [mm],      Thickness
    end
    
    methods
        function obj = Bearing(nam, typ, kx, ky, kz, ka, kb, kg, od, id, b)
            if(nargin == 0)
                nam = "-*-";      typ = "none";
                od  = 0.0;   id  = 0.0;     b  = 0.0;
                kx  = 0.0;   ky  = 0.0;     kz = 0.0;
                ka  = 0.0;   kb  = 0.0;     kg = 0.0;
            elseif(nargin < 11)
                error("Not enough input arguments.");
            end
            
            obj.name    = nam;  obj.type    = typ;  
            obj.OD      = od;   obj.ID      = id;   obj.B       = b;
            obj.K_x     = kx;   obj.K_y     = ky;   obj.K_z     = kz;
            obj.K_alpha = ka;   obj.K_beta  = kb;   obj.K_gamma = kg;
        end
        
        function tab = disp(obj)
            tab_set = {"Type",                             "-+-",     "-",       obj.type;
                       "Translational Stiffness",          "K_x",     "N/m",     obj.K_x;
                       "Translational Stiffness",          "K_y",     "N/m",     obj.K_y;
                       "Translational Stiffness",          "K_z",     "N/m",     obj.K_z;
                       "Rotational Stiffness (rot. axis)", "K_alpha", "N-m/rad", obj.K_alpha;
                       "Rotational Stiffness",             "K_beta",  "N-m/rad", obj.K_beta;
                       "Rotational Stiffness",             "K_gamma", "N-m/rad", obj.K_gamma;
                       "Outer diameter",                   "OD",      "mm",      obj.OD;
                       "Inner diameter",                   "ID",      "mm",      obj.ID;
                       "Width",                            "B",       "mm",      obj.B;
                       };
            
            Parameter = tab_set(:, 1);
            Symbol    = tab_set(:, 2);
            Unit      = tab_set(:, 3);
            Value     = tab_set(:, 4:end);
            
            n = numel(obj);
            tab_val = cell(n + 2, 1);
            
            tab_val{1} = table(Parameter, Symbol);
            
            for idx = 1:length(obj)
                tab_val{idx + 1} = table(Value(:,idx), 'variableNames', obj(idx).name);
            end
            
            tab_val{idx + 2} = table(Unit);
            
            tab = [tab_val{:}];

            if(nargout == 0)
                disp(tab);
                clear tab;
            end
        end
    end
    
    %% Calculation:
    methods
        function val = series_association(obj)
            n = numel(obj);
            if(n < 2)
                warninig("Only one bearing.");
                val = obj;
            else
                nam = "series";
                
                kx = 0.0;   ky = 0.0;   kz = 0.0;
                ka = 0.0;   kb = 0.0;   kg = 0.0;
                for idx = 1:nargin
                    nam = join([nam obj(idx).name], " / ");

                    kx = kx + 1.0/obj(idx).K_x;
                    ky = ky + 1.0/obj(idx).K_y;
                    kz = kz + 1.0/obj(idx).K_z;
                    ka = ka + 1.0/obj(idx).K_alpha;
                    kb = kb + 1.0/obj(idx).K_beta;
                    kg = kg + 1.0/obj(idx).K_gamma;
                end
                
                kx = 1.0/kx;    ky = 1.0/ky;    kz = 1.0/kz;
                ka = 1.0/ka;    kb = 1.0/kb;    kg = 1.0/kg;
                
                val = Bearing(nam, obj(1).type, kx, ky, kz, ka, kb, kg, ...
                                   obj(1).OD, ...
                                   obj(1).ID, ...
                                   obj(1).B);
            end
        end
        
        function val = parallel_association(obj)
            n = numel(obj);
            if(n < 2)
                warninig("Only one bearing.");
                val = obj;
            else
                nam = "parallel";

                kx = 0.0;   ky = 0.0;   kz = 0.0;
                ka = 0.0;   kb = 0.0;   kg = 0.0;
                
                for idx = 1:n
                    nam = join([nam obj(idx).name], " / ");

                    kx = kx + obj(idx).K_x;
                    ky = ky + obj(idx).K_y;
                    kz = kz + obj(idx).K_z;
                    ka = ka + obj(idx).K_alpha;
                    kb = kb + obj(idx).K_beta;
                    kg = kg + obj(idx).K_gamma;
                end
                
                val = Bearing(nam, obj(1).type, kx, ky, kz, ka, kb, kg, ...
                                   obj(1).OD, ...
                                   obj(1).ID, ...
                                   obj(1).B);
            end
        end
        
        function mat = stiffness_matrix(obj)
            n = numel(obj);
            mat = zeros(6, 6, n);
            
            for idx = 1:n
                mat(:, :, idx) = diag([obj(idx).K_x ...
                                       obj(idx).K_y ...
                                       obj(idx).K_z ...
                                       obj(idx).K_alpha ...
                                       obj(idx).K_beta ...
                                       obj(idx).K_gamma]);
            end
        end
        
    end
    
    methods(Static)
%     
        function val = NREL_5MW_Bearing(stage)
            switch(stage)
                case 0 % Input shaft
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    INP_A     = Bearing("INP_A"    , "CARB", 0.0   , 1.50e10, 1.50e10, 0.0, 5.0e6, 5.0e6, 1750.0, 1250.0, 375.0);
                    INP_B     = Bearing("INP_B"    , "SRB" , 4.06e8, 1.54e10, 1.54e10, 0.0, 0.0  , 0.0  , 1220.0,  750.0, 365.0);
                    
                    val = [INP_A, INP_B];
                    
                case 1
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    Sun       = Bearing;
                    PL_A      = Bearing("PL_A"     , "CRB",  9.1e4,  9.4e9,   3.2e9,   0.0, 1.4e6, 4.5e6, 600.0,  400.0,  272.0);
                    PL_B      = Bearing("PL_B"     , "CRB",  9.1e4,  9.4e9,   3.2e9,   0.0, 1.4e6, 4.5e6, 600.0,  400.0,  272.0);
                    Ring      = Bearing;
                    PLC_A     = Bearing("PLC_A"    , "SRB",  6.6e4,  1.7e9,   1.1e9,   0.0, 5.6e5, 1.3e5, 1030.0,  710.0, 315.0);
                    PLC_B     = Bearing("PLC_B"    , "CRB",  6.6e7,  1.7e9,   1.1e9,   0.0, 5.6e5, 1.3e5, 1220.0, 1000.0, 128.0);
                    
                    val = [Sun,   Sun,  ... % Sun
                           PL_A,  PL_B, ... % Planet
                           Ring,  Ring, ... % Ring
                           PLC_A, PLC_B];   % Carrier
                    
                case 2
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    Sun       = Bearing;
                    IMS_PL_A  = Bearing("IMS_PL_A" , "CRB",  9.1e4,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 520.0,  380.0,  140.0);
                    IMS_PL_B  = Bearing("IMS_PL_B" , "CRB",  9.1e4,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 520.0,  380.0,  140.0);
                    Ring      = Bearing;
                    IMS_PLC_A = Bearing("IMS_PLC_A", "CARB", 9.1e4,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 1030.0, 710.0,  236.0);
                    IMS_PLC_B = Bearing("IMS_PLC_B", "CRB" , 9.1e7,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4,  870.0, 600.0,  155.0);
                    
                    val = [Sun,       Sun,      ... % Sun
                           IMS_PL_A,  IMS_PL_B, ... % Planet
                           Ring,      Ring,     ... % Ring
                           IMS_PLC_A, IMS_PLC_B];   % Carrier
                    
                case 3
%                                       Name,        Type,   K_x,    K_y,     K_z,     K_a, K_b,   K_g,   OD,     ID,     B
                    IMS_A     = Bearing("IMS_A"    , "CRB",  0.0,    6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4, 360.0,  200.0,   98.0);
                    IMS_B     = Bearing("IMS_B"    , "TRB",  7.4e7,  5.0e8,   5.0e8,   0.0, 1.6e6, 1.8e6, 460.0,  200.0,  100.0);
                    IMS_C     = Bearing("IMS_C"    , "TRB",  7.8e7,  7.4e8,   3.3e8,   0.0, 1.1e6, 2.5e6, 460.0,  200.0,  100.0);
                    HS_A      = Bearing("HS_A"     , "CRB",  1.3e8,  8.2e8,   8.2e8,   0.0, 1.7e5, 1.0e6, 500.0,  400.0,  100.0);
                    HS_B      = Bearing("HS_B"     , "TRB",  6.7e7,  8.0e8,   1.3e8,   0.0, 1.7e5, 1.0e6, 550.0,  410.0,   86.0);
                    HS_C      = Bearing("HS_C"     , "TRB",  8.0e7,  1.0e9,   7.3e7,   0.0, 1.7e5, 1.0e6, 550.0,  410.0,   86.0);
                    
                    val = [HS_A,  HS_B,  HS_C, ... % Pinion
                           IMS_A, IMS_B, IMS_C];   % Wheel

                otherwise
                    error("prog:input", "Option [%d] is NOT valid.", stage);
            end
        end
        
    end
end


