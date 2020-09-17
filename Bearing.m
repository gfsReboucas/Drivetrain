classdef Bearing
    %BEARING
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
                                              % CARB: CARB toroidal roller bearing
                                                      % CRB: Cylindrical Roller Bearing
                                                             % SRB: Spherical Roller Bearing
                                                                    % TRB: Tapered Roller Bearing
        x;       % [mm],      bearing axial position w. r. t. shaft 
        K_x;     % [N/m],     Translational stiffness, x axis
        K_y;     % [N/m],     Translational stiffness, y axis
        K_z;     % [N/m],     Translational stiffness, z axis
        K_alpha; % [N-m/rad], Torsional stiffness, x axis (rot.)
        K_beta;  % [N-m/rad], Torsional stiffness, y axis
        K_gamma; % [N-m/rad], Torsional stiffness, z axis

        C_x;     % [N-s/m],     Translational damping, x axis
        C_y;     % [N-s/m],     Translational damping, y axis
        C_z;     % [N-s/m],     Translational damping, z axis
        C_alpha; % [N-m-s/rad], Torsional damping, x axis (rot.)
        C_beta;  % [N-m-s/rad], Torsional damping, y axis
        C_gamma; % [N-m-s/rad], Torsional damping, z axis

        OD;      % [mm],      Outer diameter
        ID;      % [mm],      Inner diameter
        B;       % [mm],      Thickness
    end
    
    methods
%         function obj = Bearing(nam, typ, kx, ky, kz, ka, kb, kg, od, id, b)
        function obj = Bearing(varargin)%nam, typ, kx, ky, kz, ka, kb, kg, od, id, b)
            default = {'name'  , 'no_name', ...
                      'type'   , 'none', ...
                      'x'      , 0.0, ...
                      'OD'     , 0.0, ...
                      'ID'     , 0.0, ...
                      'B'      , 0.0, ...
                      'K_x'    , 0.0, ...
                      'K_y'    , 0.0, ...
                      'K_z'    , 0.0, ...
                      'K_alpha', 0.0, ...
                      'K_beta' , 0.0, ...
                      'K_gamma', 0.0, ...
                      'C_x'    , 0.0, ...
                      'C_y'    , 0.0, ...
                      'C_z'    , 0.0, ...
                      'C_alpha', 0.0, ...
                      'C_beta' , 0.0, ...
                      'C_gamma', 0.0};
                  
            default = scaling_factor.process_varargin(default, varargin);
            
            obj.name    = default.name;
            obj.type    = default.type;
            obj.OD      = default.OD;
            obj.ID      = default.ID;
            obj.B       = default.B;
            
            obj.x       = default.x;
            obj.K_x     = default.K_x;
            obj.K_y     = default.K_y;
            obj.K_z     = default.K_z;
            obj.K_alpha = default.K_alpha;
            obj.K_beta  = default.K_beta;
            obj.K_gamma = default.K_gamma;
            
            obj.C_x     = default.C_x;
            obj.C_y     = default.C_y;
            obj.C_z     = default.C_z;
            obj.C_alpha = default.C_alpha;
            obj.C_beta  = default.C_beta;
            obj.C_gamma = default.C_gamma;
        end
        
        function tab = disp(obj)
            if(isempty(obj))
                disp("\t0x0 empty Bearing object")
            else
                tab_set = {"Type",                             "-+-",     "-",         obj.type;
                           "Translational Stiffness",          "K_x",     "N/m",       obj.K_x;
                           "Translational Stiffness",          "K_y",     "N/m",       obj.K_y;
                           "Translational Stiffness",          "K_z",     "N/m",       obj.K_z;
                           "Rotational Stiffness (rot. axis)", "K_alpha", "N-m/rad",   obj.K_alpha;
                           "Rotational Stiffness",             "K_beta",  "N-m/rad",   obj.K_beta;
                           "Rotational Stiffness",             "K_gamma", "N-m/rad",   obj.K_gamma;
                           "Translational Damping",            "C_x",     "N-s/m",     obj.C_x;
                           "Translational Damping",            "C_y",     "N-s/m",     obj.C_y;
                           "Translational Damping",            "C_z",     "N-s/m",     obj.C_z;
                           "Rotational Damping (rot. axis)",   "C_alpha", "N-m-s/rad", obj.C_alpha;
                           "Rotational Damping",               "C_beta",  "N-m-s/rad", obj.C_beta;
                           "Rotational Damping",               "C_gamma", "N-m-s/rad", obj.C_gamma;
                           "Outer diameter",                   "OD",      "mm",        obj.OD;
                           "Inner diameter",                   "ID",      "mm",        obj.ID;
                           "Width",                            "B",       "mm",        obj.B;
                           };

                Parameter = tab_set(:, 1);
                Symbol    = tab_set(:, 2);
                Unit      = tab_set(:, 3);
                Value     = tab_set(:, 4:end);
                
                n = numel(obj);
                
                tab_val = cell(n + 2, 1);
                
                tab_val{1} = table(Parameter, Symbol);
                
                for idx = 1:n
                    tab_val{idx + 1} = table(Value(:, idx), 'variableNames', obj(idx).name);
                end
                
                tab_val{idx + 2} = table(Unit);
                
                tab = [tab_val{:}];
                
                if(nargout == 0)
                    disp(tab);
                    clear tab;
                end
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

                cx = 0.0;   cy = 0.0;   cz = 0.0;
                ca = 0.0;   cb = 0.0;   cg = 0.0;

                for idx = 1:nargin
                    nam = join([nam obj(idx).name], " / ");

                    kx = kx + 1.0/obj(idx).K_x;
                    ky = ky + 1.0/obj(idx).K_y;
                    kz = kz + 1.0/obj(idx).K_z;
                    ka = ka + 1.0/obj(idx).K_alpha;
                    kb = kb + 1.0/obj(idx).K_beta;
                    kg = kg + 1.0/obj(idx).K_gamma;
                    
                    cx = cx + 1.0/obj(idx).C_x;
                    cy = cy + 1.0/obj(idx).C_y;
                    cz = cz + 1.0/obj(idx).C_z;
                    ca = ca + 1.0/obj(idx).C_alpha;
                    cb = cb + 1.0/obj(idx).C_beta;
                    cg = cg + 1.0/obj(idx).C_gamma;
                end
                
                if(1.0/kx == inf),  kx = 0.0;  else,  kx = 1.0/kx;  end
                if(1.0/ky == inf),  ky = 0.0;  else,  ky = 1.0/ky;  end
                if(1.0/kz == inf),  kz = 0.0;  else,  kz = 1.0/kz;  end
                if(1.0/ka == inf),  ka = 0.0;  else,  ka = 1.0/ka;  end
                if(1.0/kb == inf),  kb = 0.0;  else,  kb = 1.0/kb;  end
                if(1.0/kg == inf),  kg = 0.0;  else,  kg = 1.0/kg;  end

                if(1.0/cx == inf),  cx = 0.0;  else,  cx = 1.0/cx;  end
                if(1.0/cy == inf),  cy = 0.0;  else,  cy = 1.0/cy;  end
                if(1.0/cz == inf),  cz = 0.0;  else,  cz = 1.0/cz;  end
                if(1.0/ca == inf),  ca = 0.0;  else,  ca = 1.0/ca;  end
                if(1.0/cb == inf),  cb = 0.0;  else,  cb = 1.0/cb;  end
                if(1.0/cg == inf),  cg = 0.0;  else,  cg = 1.0/cg;  end
                
                val = Bearing('name'   , nam      , 'type'  , obj(1).type, ...
                              'K_x'    , kx       , 'K_y'   , ky         , 'K_z'    , kz,...
                              'K_alpha', ka       , 'K_beta', kb         , 'K_gamma', kg, ...
                              'C_x'    , cx       , 'C_y'   , cy         , 'C_z'    , cz,...
                              'C_alpha', ca       , 'C_beta', cb         , 'C_gamma', cg, ...
                              'OD'     , obj(1).OD, 'ID'    , obj(1).ID  , 'B'      , obj(1).B);
            end
        end
        
        function val = parallel_association(obj)
            n = numel(obj);
            if(n < 2)
                warning("Bearing:parallel_assoc", "Only one bearing.");
                val = obj;
            else
                nam = "parallel";

                kx = 0.0;   ky = 0.0;   kz = 0.0;
                ka = 0.0;   kb = 0.0;   kg = 0.0;

                cx = 0.0;   cy = 0.0;   cz = 0.0;
                ca = 0.0;   cb = 0.0;   cg = 0.0;

                for idx = 1:n
                    nam = join([nam obj(idx).name], " / ");

                    kx = kx + obj(idx).K_x;
                    ky = ky + obj(idx).K_y;
                    kz = kz + obj(idx).K_z;
                    ka = ka + obj(idx).K_alpha;
                    kb = kb + obj(idx).K_beta;
                    kg = kg + obj(idx).K_gamma;
                    
                    cx = cx + obj(idx).C_x;
                    cy = cy + obj(idx).C_y;
                    cz = cz + obj(idx).C_z;
                    ca = ca + obj(idx).C_alpha;
                    cb = cb + obj(idx).C_beta;
                    cg = cg + obj(idx).C_gamma;
                end
                
                val = Bearing('name'   , nam      , 'type'  , obj(1).type, ...
                              'K_x'    , kx       , 'K_y'   , ky         , 'K_z'    , kz,...
                              'K_alpha', ka       , 'K_beta', kb         , 'K_gamma', kg, ...
                              'C_x'    , cx       , 'C_y'   , cy         , 'C_z'    , cz,...
                              'C_alpha', ca       , 'C_beta', cb         , 'C_gamma', cg, ...
                              'OD'     , obj(1).OD, 'ID'    , obj(1).ID  , 'B'      , obj(1).B);
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
        
        function mat = damping_matrix(obj)
            n = numel(obj);
            mat = zeros(6, 6, n);
            
            for idx = 1:n
                mat(:, :, idx) = diag([obj(idx).C_x ...
                                       obj(idx).C_y ...
                                       obj(idx).C_z ...
                                       obj(idx).C_alpha ...
                                       obj(idx).C_beta ...
                                       obj(idx).C_gamma]);
            end
        end
        
    end
    
end


