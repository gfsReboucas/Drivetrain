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
            if(isempty(obj))
                disp("\t0x0 empty Bearing object")
            else
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
    
end


