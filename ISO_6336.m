classdef ISO_6336 < Gear_Set
    %ISO_6336 do some calculations following the ISO 6336 standard.
    %
    % References:
    %   [1] ISO 6336-1:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 1: Basic principles, introduction and general 
    % influence factors.
    %   [2] ISO 6336-2:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 2: Calculation of surface durability (pitting).
    %   [3] ISO/TR 6336-30:2017 Calculation of load capacity of spur and 
    % helical gears -- Calculation examples for the application of ISO 6336 
    % parts 1, 2, 3, 5
    %   [4] Arnaudov, K., Karaivanov, D. (2019). Planetary Gear Trains. 
    % Boca Raton: CRC Press, https://doi.org/10.1201/9780429458521
    %
    
    properties(Access = public)
        P;                        % [kW], Rated power
        n_out;                    % [1/min.], Output speed
%         gear_set (1, :) Gear_Set; % [-], Gear stage
        S_Hmin;                   % [-], Minimum required safety factor for surface durability (pitting)
        S_Fmin;                   % [-], Minimum required safety factor for tooth bending strength
        L_h;                      % [h], Required life
        K_A;                      % [-], Application factor
        calculation;              % [-], calculation method: by [KISSsoft] or [default]
        S_H;                      % [-], Pitting safety factor
        S_F;                      % [-], Tooth bending safety factor
    end
    
    properties(Dependent)
        K_gamma; % [-], Mesh load factor
    end
    
    methods
        function obj = ISO_6336(gset, calc, PP, nout, SHm, SFm, Lh, KA)
            
            if(~strcmp(class(gset), "Gear_Set"))
                error('gset is not a Gear_Set object.');
            end
            
            obj@Gear_Set(gset.configuration, ...
                         gset.m_n, ...
                         gset.alpha_n, ...
                         gset.z, ...
                         gset.b, ...
                         gset.x, ...
                         gset.beta, ...
                         gset.k, ...
                         gset.bore_ratio, ...
                         gset.N_p, ...
                         gset.a_w, ...
                         gset.type, ...
                         gset.bearing, ...
                         gset.out_shaft, ...
                         gset.Q);

            obj.calculation = calc;
            
            if(strcmpi(obj.calculation, 'KISSsoft') && ~KISSsoftCOM.is_installed())
                obj.calculation = 'default';
                warning('KISSsoft COM interface not found. Using [default] calculation method.');
            end
            
            obj.P        = PP;
            obj.n_out    = nout;
            obj.S_Hmin   = SHm;
            obj.S_Fmin   = SFm;
            obj.L_h      = Lh;
            obj.K_A      = KA;
        end
    end
    %% Calculation methods:
    methods
        function [SH, SF] = safety_factors(obj)
            if(strcmpi(obj.calculation, 'KISSsoft'))
                [SH, SF] = obj.run_KS();
            elseif(strcmpi(obj.calculation, 'default'))
                warning('Calculation of the tooth bending safety factor not implemented yet.');
                SH = obj.Pitting();
                SF = zeros(size(SH));
            end
        end
        
        function [SH, SF] = run_KS(obj)
            ks = actxserver('KISSsoftCOM.KISSsoft');
            ks.SetSilentMode(true);
            
            if(strcmp(obj.configuration, "parallel"))
                ks.GetModule('Z012', false);
                std_file = 'CylGearPair 1 (spur gear).Z12';
                geo_meth = false;
            elseif(strcmp(obj.configuration, "planetary"))
                ks.GetModule('Z014', false);
                std_file = 'PlanetarySet 1 (ISO6336).Z14';
                geo_meth = true;
            end
            
            file_name = sprintf('C:\\Program Files (x86)\\KISSsoft 03-2017\\example\\%s', std_file);
            
            try
                ks.LoadFile(file_name);
            catch err
                ks.ReleaseModule();
                error(err.identifier, "%s", err.message);
                return;
            end
            
            ks.SetVar('ZS.AnzahlZwi', num2str(        obj.N_p));              % number of planets
            ks.SetVar('ZS.Geo.mn'   , num2str(        obj.m_n     , '%.6f')); % normal module
            ks.SetVar('ZP[0].a'     , num2str(        obj.a_w     , '%.6f')); % center distance
            ks.SetVar('ZS.Geo.alfn' , num2str(deg2rad(obj.alpha_n), '%.6f')); % normal pressure angle
            ks.SetVar('ZS.Geo.beta' , num2str(deg2rad(obj.beta)   , '%.6f')); % helix angle
            
            ks.SetVar('RechSt.GeometrieMeth', num2str(geo_meth));    % tooth geometry according to ISO 21771:2007
            
        end
        
        function SH = Pitting(obj, R_a, R_z, nu_40, line, C_a)
            % preparatory calculations:
            T_1 = (obj.P*1.0e3)/(obj.n_out*pi/30.0);
            
            T_1 = abs(T_1);
            
            SH = obj.calculate_SH(T_1, C_a, E, nu, rho, ...
                                  sigma_Hlim, line, nu_40, R_z);
        end
        
        function SH = calculate_SH(obj, T_1)
        end
    end
    
    %% Set methods:
    
    %% Get methods:
    methods
        function val = get.K_gamma(obj)
            Np = obj.gear_set.N_p;
            if(Np == 3)
                val = 1.1;
            elseif(Np == 4)
                val = 1.25;
            elseif(Np == 5)
                val = 1.35;
            elseif(Np == 6)
                val = 1.44;
            elseif(Np == 7)
                val = 1.47;
            else
                val = 1.0;
            end
        end
    end
end
