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
            end
            
            obj.name    = nam;  obj.type    = typ;  
            obj.OD      = od;   obj.ID      = id;   obj.B       = b;
            obj.K_x     = kx;   obj.K_y     = ky;   obj.K_z     = kz;
            obj.K_alpha = ka;   obj.K_beta  = kb;   obj.K_gamma = kg;
        end
        
        function tab = disp(obj)
            tab_set = {"Name",                             "-+-",     "-",       obj.name;
                       "Type",                             "-+-",     "-",       obj.type;
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
            
            Parameter = tab_set(:,1);
            Symbol    = tab_set(:,2);
            Unit      = tab_set(:,3);
            Value     = tab_set(:,4);

            tab = table(Parameter, Symbol, Value, Unit, ...
                'variableNames', ["Parameter", "Symbol", "Value", "Unit"]);
            
            if(nargout == 0)
                disp(tab);
                clear tab;
            end
        end
    end
    
    %% Calculation:
    methods
        function obj = series_association(varargin)
            if(nargin < 2)
                error("prog:input", "Minimum two bearings");
            end
            
            nam = "series";
            kx = 0.0;   ky = 0.0;   kz = 0.0;
            ka = 0.0;   kb = 0.0;   kg = 0.0;
            for idx = 1:nargin
                nam = join([nam varargin{idx}.name], " / ");
                
                kx = kx + 1.0/varargin{idx}.K_x;
                ky = ky + 1.0/varargin{idx}.K_y;
                kz = kz + 1.0/varargin{idx}.K_z;
                ka = ka + 1.0/varargin{idx}.K_alpha;
                kb = kb + 1.0/varargin{idx}.K_beta;
                kg = kg + 1.0/varargin{idx}.K_gamma;
            end
            
            kx = 1.0/kx;    ky = 1.0/ky;    kz = 1.0/kz;
            ka = 1.0/ka;    kb = 1.0/kb;    kg = 1.0/kg;

            obj = Bearing(nam, varargin{1}.type, kx, ky, kz, ka, kb, kg, ...
                               varargin{1}.OD, ...
                               varargin{1}.ID, ...
                               varargin{1}.B);
        end
        
        function obj = parallel_association(varargin)
            if(nargin < 2)
                error("prog:input", "Minimum two bearings");
            end
            
            nam = "parallel";
            
            kx = 0.0;   ky = 0.0;   kz = 0.0;
            ka = 0.0;   kb = 0.0;   kg = 0.0;
            
            for idx = 1:nargin
                nam = join([nam varargin{idx}.name], " / ");
                
                kx = kx + varargin{idx}.K_x;
                ky = ky + varargin{idx}.K_y;
                kz = kz + varargin{idx}.K_z;
                ka = ka + varargin{idx}.K_alpha;
                kb = kb + varargin{idx}.K_beta;
                kg = kg + varargin{idx}.K_gamma;
            end
            
            obj = Bearing(nam, varargin{1}.type, kx, ky, kz, ka, kb, kg, ...
                               varargin{1}.OD, ...
                               varargin{1}.ID, ...
                               varargin{1}.B);
        end
        
        function mat = stiffness_matrix(obj)
            mat = diag([obj.K_x ...
                        obj.K_y ...
                        obj.K_z ...
                        obj.K_alpha ...
                        obj.K_beta ...
                        obj.K_gamma]);
        end
        
    end
    
    methods(Static)
%     
%         function val = NREL_5MW_Bearing(stage)
%             switch(stage)
%                 case 0 % Input shaft
%                     INP_A     = Bearing("INP_A",     0.0,    1.5e10,  1.5e10,  0.0, 5.0e6, 5.0e6);
%                     INP_B     = Bearing("INP_B",     4.06e8, 1.54e10, 1.54e10, 0.0, 0.0,   0.0  );
%                     
%                     val = [INP_A, INP_B];
%                     
%                 case 1
%                     Sun       = Bearing;
%                     PLC_A     = Bearing("PLC_A",     6.6e4,  1.7e9,   1.1e9,   0.0, 5.6e5, 1.3e5);
%                     PLC_B     = Bearing("PLC_B",     6.6e7,  1.7e9,   1.1e9,   0.0, 5.6e5, 1.3e5);
%                     PL_A      = Bearing("PL_A",      9.1e4,  9.4e9,   3.2e9,   0.0, 1.4e6, 4.5e6);
%                     PL_B      = Bearing("PL_B",      9.1e4,  9.4e9,   3.2e9,   0.0, 1.4e6, 4.5e6);
%                     Ring      = Bearing;
%                     
%                     val = [Sun,   Sun;    % Sun
%                            PL_A,  PL_B;   % Planet
%                            Ring,  Ring;   % Ring
%                            PLC_A, PLC_B]; % Carrier
%                     
%                 case 2
%                     Sun       = Bearing;
%                     IMS_PLC_A = Bearing("IMS_PLC_A", 9.1e4,  6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4);
%                     IMS_PLC_B = Bearing("IMS_PLC_B", 9.1e7, 6.0e7,    1.2e9,   0.0, 7.5e4, 7.5e4);
%                     IMS_PL_A  = Bearing("IMS_PL_A",  9.1e4, 6.0e7,    1.2e9,   0.0, 7.5e4, 7.5e4);
%                     IMS_PL_B  = Bearing("IMS_PL_B",  9.1e4, 6.0e7,    1.2e9,   0.0, 7.5e4, 7.5e4);
%                     Ring      = Bearing;
%                     
%                     val = [Sun,       Sun;        % Sun
%                            IMS_PL_A,  IMS_PL_B;   % Planet
%                            Ring,      Ring;       % Ring
%                            IMS_PLC_A, IMS_PLC_B]; % Carrier
%                     
%                 case 3
%                     IMS_A     = Bearing("IMS_A",     0.0,    6.0e7,   1.2e9,   0.0, 7.5e4, 7.5e4);
%                     IMS_B     = Bearing("IMS_B",     7.4e7,  5.0e8,   5.0e8,   0.0, 1.6e6, 1.8e6);
%                     IMS_C     = Bearing("IMS_C",     7.8e7,  7.4e8,   3.3e8,   0.0, 1.1e6, 2.5e6);
%                     HS_A      = Bearing("HS_A",      1.3e8,  8.2e8,   8.2e8,   0.0, 1.7e5, 1.0e6);
%                     HS_B      = Bearing("HS_B",      6.7e7,  8.0e8,   1.3e8,   0.0, 1.7e5, 1.0e6);
%                     HS_C      = Bearing("HS_C",      8.0e7,  1.0e9,   7.3e7,   0.0, 1.7e5, 1.0e6);
%                     
%                     val = [IMS_A, IMS_B, IMS_C; % Pinion
%                            HS_A,  HS_B,  HS_C]; % Wheel
% 
%                 otherwise
%                     error("prog:input", "Option [%d] is NOT valid.", stage);
%             end
%         end
%         
    end
end


