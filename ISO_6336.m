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
    
    properties(Access = private)
        gear_set (1, :) Gear_Set; % [-], Gear stage
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

            obj.gear_set = gset;

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
                [SH, SF] = obj.KISSsoft_calculation();
            elseif(strcmpi(obj.calculation, 'default'))
                warning('Calculation of the tooth bending safety factor not implemented yet.');
                SH = obj.Pitting();
                SF = zeros(size(SH));
            end
        end
        
        function [SH, SF] = KISSsoft_calculation(obj, varargin)
            if(nargin <= 2)
                lubricant_ID = '11220';
                stage_idx    = 0;
                save_report  = false;
                show_report  = false;
            else
                lubricant_ID = varargin{1};
                stage_idx    = varargin{2};
                save_report  = varargin{3}; % true
                show_report  = varargin{4}; % true
            end
            
            ks = obj.KISSsoft();
            ks.set_var('ZS.P'       , obj.P); % rated power
            ks.set_var('ZS.Pnominal', obj.P); % rated power
            
            ks.set_var('ZS.Oil.SchmierTypID', lubricant_ID);
            
            flag_impul = true;
            if(strcmp(obj.configuration, "parallel"))
                flag_impul = true;
                template_code = 12;
            elseif(strcmp(obj.configuration, 'planetary'))
                n_c = obj.n_out/obj.u;
                T_c = obj.P/(n_c*pi/30.0);
                n_planet = n_c*(1.0 - obj.z(3)./obj.z(2));
                ks.set_var('ZR[1].n'               , abs(n_planet)); % speed at gear 2: planet
                ks.set_var('ZR[2].n'               , 0.0          ); % speed at gear 3: ring
                ks.set_var('ZS.Planet.nStegnominal', n_c          ); % speed at planet carrier
                ks.set_var('ZS.Planet.nSteg'       , n_c          ); % speed at planet carrier
                
                ks.set_var('ZR[1].T'               , abs(n_planet)); % speed at gear 2: planet
                ks.set_var('ZS.Planet.TStegnominal', T_c          ); % torque at planet carrier
                ks.set_var('ZS.Planet.nSteg'       , T_c          ); % torque at planet carrier
                
                flag_impul = false;
                template_code = 14;
            end
            
            n_1 = obj.n_out;
            T_1 = obj.P/(n_1*pi/30.0);
            ks.set_var('ZR[0].Tnominal', T_1       ); % torque at gear 1
            ks.set_var('ZR[0].n'       , n_1       ); % speed at gear 1
            ks.set_var('ZP[0].Impuls'  , flag_impul); % sense of rotation, gear 1

            ks.set_var('ZS.SSi.Flanke'      , obj.S_Hmin);
            ks.set_var('ZS.SSi_0to05.Flanke', obj.S_Hmin);
            ks.set_var('ZS.SSi_at1.Flanke'  , obj.S_Hmin);
            ks.set_var('ZS.SSi_2to20.Flanke', obj.S_Hmin);
            ks.set_var('ZS.SSi_fix.Flanke'  , obj.S_Hmin);
            
            ks.set_var('ZS.SSi.Fuss'      , obj.S_Fmin);
            ks.set_var('ZS.SSi_0to05.Fuss', obj.S_Fmin);
            ks.set_var('ZS.SSi_at1.Fuss'  , obj.S_Fmin);
            ks.set_var('ZS.SSi_2to20.Fuss', obj.S_Fmin);
            ks.set_var('ZS.SSi_fix.Fuss'  , obj.S_Fmin);
            
            ks.set_var('ZS.H'   , obj.L_h);
            ks.set_var('ZS.KA'  , obj.K_A);
            ks.set_var('ZS.Kgam', obj.K_gamma);
            
            if(~ks.calculate())
                delete(ks);
            end
            
            kv = ks.get_var('ZS.KVcalc');
            
            flag_kv = false;
            if(kv < 1.05)
                flag_kv = true;
                for idx = 1:length(obj.z)
                    ks.set_var(sprintf('ZP[%d].KV.KV', idx), 1.05);
                end
                
                ks.set_var('Zst.KVFlag', flag_kv);
            end
            
            flag_khb = false;
            for idx = 1:2
                khb = ks.get_var(sprintf('ZP[%d].KHb', idx - 1));
                
                if(khb < 1.15)
                    flag_khb = true;
                    ks.set_var(sprintf('ZP[%d].KHb', idx - 1), 1.15);
                end
            end
            
            if(flag_khb)
                ks.set_var('Zst.KHbVariant', true);
            end
            
            if(~ks.calculate())
                delete(ks);
            end
            
            if(flag_kv || flag_khb)
                ks.calculate();
            end
            
            max_Np = obj.N_p + 1;
            
            SH = zeros(max_Np, 1);
            SF = zeros(max_Np, 1);
            
            for idx = 1:length(obj.z)
                SH(idx) = ks.get_var(sprintf('ZPP[%d].Flanke.SH', idx - 1));
                SF(idx) = ks.get_var(sprintf('ZPP[%d].Fuss.SF'  , idx - 1));
            end
            
            SH(SH == 0) = [];
            SF(SF == 0) = [];
            
            if(isrow(SH))
                SH = SH';
            end
            
            if(isrow(SF))
                SF = SF';
            end
            
            if(save_report)
                if(strcmpi(class(obj), "Drivetrain"))
                    report_name = sprintf("%s\\stage_%02d.rtf" , pwd, stage_idx);
                    file_name   = sprintf("%s\\stage_%02d.Z0%d", pwd, stage_idx, template_code);
                else
                    report_name = sprintf("%s\\@%s\\stage_%02d.rtf" , pwd, class(obj), stage_idx);
                    file_name   = sprintf("%s\\@%s\\stage_%02d.Z0%d", pwd, class(obj), stage_idx, template_code);
                end
                
                ks.save_file(file_name);
                version = ks.get_version;
                template = sprintf("C:\\Program Files (x86)\\KISSsoft %s\\rpt\\Z0%dLe0.RPT", version, template_code);
                
                ks.write_report(template, ... % template file
                    report_name, show_report, 0); % output format
            else
                ks.report(show_report);
            end
            
            delete(ks);
        end
        
        function SH = Pitting(obj, varargin)
            
            if(nargin == 1)
                nu_40 = 220.0;
                line  = 4;
                C_a   = 0.0;
            else
                nu_40 = varargin{1};
                line  = varargin{2};
                C_a   = varargin{3};
            end
            
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
