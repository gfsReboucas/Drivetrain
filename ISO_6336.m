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
        P_rated;     % [kW], Rated power
        n_out;       % [1/min.], Output speed
        S_Hmin;      % [-], Minimum required safety factor for surface durability (pitting)
        S_Fmin;      % [-], Minimum required safety factor for tooth bending strength
        L_h;         % [h], Required life
        K_A;         % [-], Application factor
        calculation; % [-], calculation method: [KISSsoft] or [default]
        S_H;         % [-], Pitting safety factor
        S_F;         % [-], Tooth bending safety factor
    end
    
    properties(Dependent)
        K_gamma;      % [-], Mesh load factor
        T_1;          % [N-m], Applied torque
        sigma_HP_ref; % [N/mm^2], Permissive contact stress (reference)
    end
    
    properties(Access = private)
        gear_set (1, :) Gear_Set; % [-], Gear stage
    end
    
    methods
        function obj = ISO_6336(gset, varargin)
            if(~exist('gset', "var"))
                gset = Gear_Set();
            elseif(~strcmp(class(gset), "Gear_Set"))
                error('gset is not a Gear_Set object.');
            end
            
            T_1 = 9.0e3;
            n_1 = 360.0;
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            default = {'calculation', 'default', ...
                       'P_rated'    , P_r, ...
                       'n_out'      , n_1, ...
                       'S_Hmin'     , 1.0, ...
                       'S_Fmin'     , 1.0, ...
                       'L_h'        , 50.0e3, ...
                       'K_A'        , 1.0};
            
            default = process_varargin(default, varargin);
            
            obj@Gear_Set('configuration', gset.configuration, ...
                         'm_n'          , gset.m_n, ...
                         'alpha_n'      , gset.alpha_n, ...
                         'z'            , gset.z, ...
                         'b'            , gset.b, ...
                         'x'            , gset.x, ...
                         'beta'         , gset.beta, ...
                         'k'            , gset.k, ...
                         'bore_ratio'   , gset.bore_ratio, ...
                         'N_p'          , gset.N_p, ...
                         'a_w'          , gset.a_w, ...
                         'rack_type'    , gset.type, ...
                         'bearing'      , gset.bearing, ...
                         'shaft'        , gset.output_shaft, ...
                         'Q'            , gset.Q, ...
                         'R_a'          , gset.R_a);

            obj.gear_set = gset;

            if(strcmpi(default.calculation, 'KISSsoft'))
                obj.calculation = 'KISSsoft';
                if(~KISSsoftCOM.is_installed())
                    obj.calculation = 'default';
                    warning('ISO_6336:KS', ['KISSsoft\x00a9 COM interface not found. ', ...
                                            'Using [DEFAULT] calculation method.']);
                end
            elseif(strcmpi(default.calculation, 'default'))
                obj.calculation = 'default';
            end
            
            obj.P_rated  = default.P_rated;
            obj.n_out    = default.n_out;
            obj.S_Hmin   = default.S_Hmin;
            obj.S_Fmin   = default.S_Fmin;
            obj.L_h      = default.L_h;
            obj.K_A      = default.K_A;
        end
        
        function disp(obj)
            % to be done...
        end
    end
    
    %% Calculation methods:
    methods
        function [SH, SF] = safety_factors(obj, varargin)
            if(strcmpi(obj.calculation, 'KISSsoft'))
                [SH, SF] = obj.KISSsoft_calculation(varargin{:});
            elseif(strcmpi(obj.calculation, 'default'))
                warning('ISO_6336:SF', 'Calculation of the tooth bending safety factor not implemented yet.');
                SH = obj.Pitting(varargin{:});
                SF = nan(size(SH));
            end
        end
        
        function [SH, SF] = KISSsoft_calculation(obj, varargin)
            default = {'lubricant_ID', '11220', ...
                       'stage_idx'   ,  0, ...
                       'save_report' ,  false, ...
                       'show_report' ,  false};
            
            default = process_varargin(default, varargin);
            
            ks = obj.KISSsoft();
            ks.set_var('ZS.P'       , obj.P_rated); % rated power
            ks.set_var('ZS.Pnominal', obj.P_rated); % rated power
            
            ks.set_var('ZS.Oil.SchmierTypID', default.lubricant_ID);
            
            flag_impul = true;
            if(strcmp(obj.configuration, "parallel"))
                flag_impul = true;
                template_code = 12;
            elseif(strcmp(obj.configuration, 'planetary'))
                n_c = obj.n_out/obj.u;
                T_c = obj.P_rated/(n_c.*pi/30.0);
                n_planet = n_c.*(1.0 - obj.z(3)./obj.z(2));
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
            
            ks.set_var('ZR[0].Tnominal', obj.T_1   ); % torque at gear 1
            ks.set_var('ZR[0].n'       , obj.n_out ); % speed at gear 1
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
                for idx = 1:2
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
            
            SH(SH == 0.0) = [];
            SF(SF == 0.0) = [];
            
            SH(isnan(SH)) = [];
            SF(isnan(SF)) = [];
            
            if(isrow(SH))
                SH = SH';
            end
            
            if(isrow(SF))
                SF = SF';
            end
            
            if(default.save_report)
                if(strcmpi(class(obj), "Drivetrain"))
                    report_name = sprintf("%s\\stage_%02d.rtf" , pwd, default.stage_idx);
                    file_name   = sprintf("%s\\stage_%02d.Z0%d", pwd, default.stage_idx, template_code);
                else
                    report_name = sprintf("%s\\@%s\\stage_%02d.rtf" , pwd, class(obj), default.stage_idx);
                    file_name   = sprintf("%s\\@%s\\stage_%02d.Z0%d", pwd, class(obj), default.stage_idx, template_code);
                end
                
                ks.save_file(file_name);
                version = ks.get_version;
                template = sprintf("C:\\Program Files (x86)\\KISSsoft %s\\rpt\\Z0%dLe0.RPT", version, template_code);
                
                ks.write_report(template, ... % template file
                    report_name, default.show_report, 0); % output format
            else
                ks.report(default.show_report);
            end
        end
        
        function [SH, out] = Pitting(obj, varargin)
            % [N/mm^2],  Allowable contact stress number:
            sigma_Hlim = Material().sigma_Hlim*1.0e-6;
            % Tip relief by running-in, Table 8, [1]:
            C_ay = (1.0/18.0)*(sigma_Hlim/97.0 - 18.45)^2 + 1.5;
            
            default = {'nu_40', 220.0, ...
                       'line' ,   4, ...
                       'C_a'  , C_ay};
            
            default = process_varargin(default, varargin);
            
            % [N/mm^2],  Young's modulus:
            E          = Material().E*1.0e-6;
            % [-],       Poisson's ratio:
            nu         = Material().nu;
            % [kg/mm^3], Density
            rho        = Material().rho*1.0e-9;
            
            out = struct();

            f_pb = max(obj.f_pt.*cosd(obj.alpha_t));
            out.f_pb = f_pb;
            
            if(f_pb >= 40.0) % [um]
                y_alpha = 3.0; % [um]
            else
                y_alpha = f_pb*75.0e-3;
            end
            out.y_alpha = y_alpha;
            
            f_falpha = max(obj.f_falpha);
            out.f_falpha = f_falpha;
            % Estimated running allowance (pitch deviation):
            y_p = y_alpha;
            out.y_p = y_p;
            % Estimated running allowance (flank deviation):
            y_f = f_falpha*75.0e-3;
            out.y_f = y_f;
            
            f_pbeff = f_pb - y_p;
            out.f_pbeff = f_pbeff;
            
            f_falphaeff = f_falpha - y_f; % 0.925*y_f
            out.f_falphaeff = f_falphaeff;
            
            % Pitch line velocity:
            v = obj.pitch_line_velocity();
            out.v = v;
            
            Z_eps = obj.contact_ratio_factor();
            out.Z_eps = Z_eps;
            
            % [N], Nominal tangential load:
            F_t = 2.0e3*(obj.T_1/obj.d(1))/obj.N_p;
            F_t = abs(F_t);
            out.F_t = F_t;
            
            line_load = F_t*obj.K_A*obj.K_gamma/mean(obj.b);
            
            K_v = obj.dynamic_factor(v, rho, line_load, default.C_a, ...
            f_pbeff, f_falphaeff);
            out.K_v = K_v;
            
            K_Hbeta = obj.face_load_factor();
            out.K_Hbeta = K_Hbeta;
            
            % Determinant tangential load in a transverse plane:
            F_tH = F_t*obj.K_gamma*obj.K_A*K_v*K_Hbeta;
            out.F_tH = F_tH;
            
            term = obj.c_gamma_alpha*(f_pb - y_alpha)/(F_tH/obj.b);
            K_Halpha = obj.transverse_load_factor(term, Z_eps);
            out.K_Halpha = K_Halpha;
            
            % Zone factor: (sec. 6.1)
            Z_H = obj.zone_factor();
            out.Z_H = Z_H;
            % Single pair tooth contact factor (sec. 6.2)
            Z_BD = obj.tooth_contact_factor();
            out.Z_BD = Z_BD;
            
            % Elasticity factor: (sec. 7)
            Z_E = sqrt(E/(2.0*pi*(1.0 - nu^2)));
            out.Z_E = Z_E;
            
            % Helix angle factor: (sec. 9)
            Z_beta = 1.0/sqrt(cosd(obj.beta));
            out.Z_beta = Z_beta;
            
            [Z_L, Z_v] = obj.lubrication_velocity_factor(sigma_Hlim, default.nu_40, v);
            out.Z_L = Z_L;
            out.Z_v = Z_v;
            
            % Roughness factor:
            Z_R = obj.roughness_factor(sigma_Hlim);
            out.Z_R = Z_R;
            
            % Work hardening factor:
            Z_W = 1.0;
            out.Z_W = Z_W;
            
            % Size factor:
            Z_X = 1.0;
            out.Z_X = Z_X;

            % Number of load cycles:
            N_L1 = obj.n_out*60.0*obj.L_h; % pinion
            N_L2 = N_L1/obj.u;             % wheel
            out.N_L = [N_L1 N_L2];
            
            % Life factor:
            % line = 2
            Z_NT1 = obj.life_factor(default.line, N_L1);
            Z_NT2 = obj.life_factor(default.line, N_L2);
            Z_NT = [Z_NT1 Z_NT2];
            out.Z_NT = Z_NT;
            
            % Contact stress:
            % Nominal contact stress at pitch point:
            num = F_t*(obj.u + 1.0);
            den = obj.d(1)*obj.b*obj.u;
            sigma_H0 = Z_H*Z_E*Z_eps*Z_beta*sqrt(num/den);
            out.sigma_H0 = sigma_H0;
            
            % nominal contact stress at pitch point:
            sigma_H = Z_BD*sigma_H0*sqrt(obj.K_gamma*obj.K_A*K_v*K_Hbeta*K_Halpha); % pinion/wheel
            out.sigma_H = sigma_H;
            
            % Permissible contact stress:
            sigma_HP = sigma_Hlim*Z_NT*Z_L*Z_v*Z_R*Z_W*Z_X/obj.S_Hmin;
            out.sigma_HP = sigma_HP;
            
            % Safety factor for surface durability (against pitting):
            SH = obj.S_Hmin*sigma_HP./sigma_H; % pinion/planet
            
            if(isrow(SH))
                SH = SH';
            end
        end
        
    end
    
    methods(Static)
        function ex01 = example01()
            
            ex01 = Gear_Set('configuration', 'parallel', ...   % configuration
                            'm_n'          , 8.0, ...          % normal module
                            'alpha_n'      , 20.0, ...         % pressure angle
                            'z'            , [17 103], ...     % number of teeth
                            'b'            , 100.0, ...        % face width
                            'x'            , [0.145 0.0], ...  % profile shift coefficient
                            'beta'         , 15.8, ...         % helix angle
                            'k'            , [1 1]*0, ...      % k
                            'bore_ratio'   , [1 1]*0.5, ...    % bore ratio
                            'N_p'          , 1, ...            % number of planets
                            'a_w'          , 500.0, ...        % center distance
                            'rack_type'    , 'D', ...          % rack type
                            'bearing'      , Bearing(), ...    %
                            'shaft'        , Shaft(), ...      %
                            'Q'            , 5.0, ...          % ISO accuracy grade
                            'R_a'          , 1.0);             % surface roughness flank
            
            T_1 = 9.0e3;
            n_1 = 360.0;
            P_r = 1.0e-3*T_1*n_1*pi/30.0;
            
            calc = ISO_6336(ex01, 'calculation', 'default', ...
                                  'P_rated'    , P_r, ...
                                  'n_out'      , n_1, ...
                                  'S_Hmin'     , 1.0, ...
                                  'S_Fmin'     , 1.0, ...
                                  'L_h'        , 50.0e3, ...
                                  'K_A'        , 1.0);

            [SH, SF] = calc.safety_factors('C_a'  ,  70.0, ...
                                           'nu_40', 320.0, ...
                                           'line' ,   2);
            
        end
        
    end
    
    methods(Access = private)
        function v = pitch_line_velocity(obj)
            v = nan;
            if(strcmp(obj.configuration, 'parallel'))
                v = (pi*obj.n_out/60.0e3)*obj.d(1);
            elseif(strcmp(obj.configuration, 'planetary'))
                v = (pi*obj.n_out/60.0e3)*(obj.d(1) - obj.a_w/obj.u);
            end
        end
        
        function Z_eps = contact_ratio_factor(obj)
            eps_a = obj.eps_alpha;
            eps_b = obj.eps_beta;
            
            if(obj.beta == 0.0)
                Z_eps = sqrt((4.0 - eps_a)/3.0);
            else
                if(eps_b < 1.0)
                    Z_eps = sqrt((1.0 - eps_b)*(4.0 - eps_a)/3.0 + eps_b/eps_a);
                else
                    Z_eps = sqrt(1.0/eps_a);
                end
            end
            
        end
        
        function K_v = dynamic_factor(obj, v, rho, line_load, C_a, f_pbeff, f_falphaeff)
            z1 = obj.z(1);
            uu = obj.u;
            cond = (v*z1/100.0)*sqrt((uu.^2)/(1.0 + uu.^2));
            if(cond < 3.0) % [m/s]
                warning('ISO_6336:KV', ['Calculating K_v using method B', ...
                    ' outside of its useful range. ', ...
                'More info at the end of Sec. 6.3.2 of ISO 6336-1.']);
            end
            
            if(strcmp(obj.configuration, 'parallel'))
                % Based on Sec. 6.5.9 of ISO 6336-1, Eq. (30), assuming gears:
                % - of solid construction, and
                % - with the same density
                
                num = pi*rho*(obj.u^2)*obj.d_m(1)^4;
                den = 8.0*(obj.u^2 + 1.0)*obj.d_b(1)^2;
                m_red = num/den;
                
                % Resonance running speed, Eq. (6) [1]:
                n_E1 = 3.0e4*sqrt(obj.c_gamma_alpha/m_red)/(pi*obj.z(1)); % [1/min]
                
                % Resonance ratio, Eq. (9) [1]:
                N = obj.n_out/n_E1;
                
                K_v = obj.dynamic_factor_from_range(line_load, N, C_a, ...
                    f_pbeff, f_falphaeff);
            elseif(strcmp(obj.configuration, 'planetary'))
                m_sun = rho*(pi/8.0)*(obj.d_m(1)^4)/(obj.d_b(1)^2); % sun
                m_pla = rho*(pi/8.0)*(obj.d_m(2)^4)/(obj.d_b(2)^2); % planet
                m_red1 = (m_pla*m_sun)/(m_sun + m_pla*obj.N_p);
                m_red2 =  m_pla;
                
                % Resonance running speed:
                n_E11 = 3.0e4*sqrt(obj.cprime/m_red1)*1.0/(pi*obj.z(1)); % [1/min]
                n_E12 = 3.0e4*sqrt(obj.cprime/m_red2)*1.0/(pi*obj.z(2)); % [1/min]
                
                % Resonance ratio:
                N_1 =  obj.n_out/n_E11;
                N_2 = (obj.n_out/obj.u)/n_E12;

                K_v1 = obj.dynamic_factor_from_range(line_load, N_1, C_a, ...
                    f_pbeff, f_falphaeff);
                
                K_v2 = obj.dynamic_factor_from_range(line_load, N_2, C_a, ...
                    f_pbeff, f_falphaeff);
                
                K_v = max(K_v1, K_v2);
            end
        end
        
        function K_v = dynamic_factor_from_range(obj, line_load, N, C_a, f_pbeff, f_falphaeff)
            % Table 8, [1]:
            C_v1 = 0.32;
            C_v5 = 0.47;
            
            eps_g = obj.eps_gamma;
            if((1.0 < eps_g) && (eps_g <= 2.0))
                C_v2 = 0.34;
                C_v3 = 0.23;
                C_v4 = 0.90;
                C_v6 = 0.47;
            elseif(eps_g > 2.0)
                C_v2 =  0.57 /(eps_g - 0.3);
                C_v3 =  0.096/(eps_g - 1.56);
                C_v4 = (0.57 - 0.05*eps_g)/(eps_g - 1.44);
                C_v6 =  0.12/(eps_g - 1.74);
            end
                
            if((1.0 < eps_g) && (eps_g <= 1.5))
                C_v7 = 0.75;
            elseif((1.5 < eps_g) && (eps_g <= 2.5))
                C_v7 = 0.125*sin(pi*(eps_g - 2.0)) + 0.875;
            elseif(eps_g > 2.5)
                C_v7 = 1.0;
            end
            
            cp = obj.cprime;
            B_p = cp*f_pbeff/line_load;        % Eq. (15)
            B_f = cp*f_falphaeff/line_load;    % Eq. (16)
            B_k = abs(1.0 - cp*C_a/line_load); % Eq. (17)
            
            % Dynamic factor:
            if(line_load < 100) % [N/mm]
                N_S = 0.5 + 0.35*sqrt(line_load/100.0); % Eq. (11), [1]
            else
                N_S = 0.85; % Eq. (12) [1]
            end
            
            if(N <= N_S)
                % Section 6.5.4, Eq. (13), [1]:
                K = C_v1*B_p + C_v2*B_f + C_v3*B_k;
                K_v = N*K + 1.0; % Eq. (14)
            elseif((N_S < N) && (N <= 1.15))
                % Section 6.5.5, Eq. (20), [1]:
                K_v = C_v1*B_p + C_v2*B_f + C_v4*B_k + 1.0;
            elseif((1.15 < N) && (N < 1.5))
                K_vN115 = C_v1*B_p + C_v2*B_f + C_v4*B_k + 1.0;
                K_vN15  = C_v5*B_p + C_v6*B_f + C_v7;
                % Section 6.5.7, Eq. (22), [1]:
                K_v = K_vN15 + (K_vN115 - K_vN15)*(1.5 - N)/0.35;
            elseif(N >= 1.5)
                % Section 6.5.6, Eq. (21), [1]:
                K_v = C_v5*B_p + C_v6*B_f + C_v7;
            end
        end
        
        function [K_Hbeta, K_Fbeta] = face_load_factor(obj)
            % The face load factor (contact stress) is not calculated from
            % ISO 6336 [1], but according to [4], Eq.(17.14):
            K_Hbeta = 1.0 + 0.4*(obj.b/obj.d(1))^2;
            
            % Method B:
            % h_1 = h_aP + h_fP + k_1*m_n;
            h_1 = abs(obj.d_a(1) - obj.d_f(1))/2.0;
            bh1 = obj.b/h_1;
            h_2 = abs(obj.d_a(2) - obj.d_f(2))/2.0;
            bh2 = obj.b/h_2;
            
            bh = min(bh1, bh2);
            
            if(bh < 3.0)
                bh = 3.0;
                
                hb = 1.0/bh;
                
                N_F = 1.0/(1.0 + hb + hb.^2);
                
                % Face load factor (root stress):
                K_Fbeta = power(K_Hbeta, N_F);
            end
        end
        
        function [K_Halpha, K_Falpha] = transverse_load_factor(obj, term, Z_eps)
            % Section 8.3.2 of [1]:
            if(obj.eps_gamma <= 2.0) % Eq. (73)
                K_Falpha = (0.9 + 0.4*term)*(obj.eps_gamma/2.0);
            else % Eq. (74)
                K_Falpha = 0.9 + 0.4*sqrt(2.0*(obj.eps_gamma - 1.0)/obj.eps_gamma)*term;
            end
            
            K_Halpha = K_Falpha;
            
            % Section 8.3.4 of [1]:
            % Transverse load factor (contact stress), Eq. (75):
            K_Halpha_lim = obj.eps_gamma/(obj.eps_alpha*Z_eps.^2);
            if(K_Halpha > K_Halpha_lim)
                K_Halpha = K_Halpha_lim;
            elseif(K_Halpha < 1.0)
                K_Halpha = 1.0;
            end
            
            % Section 8.3.5 of [1]:
            % Transverse load factor (root stress), Eq. (76):
            K_Falpha_lim = obj.eps_gamma/(0.25*obj.eps_alpha + 0.75);
            if(K_Falpha > K_Falpha_lim)
                K_Falpha = K_Falpha_lim;
            elseif(K_Falpha < 1.0)
                K_Falpha = 1.0;
            end
            
        end
        
        function Z_H = zone_factor(obj)
            num = 2.0*cosd(obj.beta_b)*cosd(obj.alpha_wt);
            den = sind(obj.alpha_wt)*cosd(obj.alpha_t).^2;
            Z_H = sqrt(num/den);
        end
        
        function Z_BD = tooth_contact_factor(obj)
            M_1 = tand(obj.alpha_wt)/sqrt((sqrt((obj.d_a(1)/obj.d_b(1))^2 - 1.0) - ...
                2.0*pi/obj.z(1))*(sqrt((obj.d_a(2)/obj.d_b(2))^2 - 1.0) - ...
                (obj.eps_alpha - 1.0)*2.0*pi/obj.z(2)));
            M_2 = tan(obj.alpha_wt)/sqrt((sqrt((obj.d_a(2)/obj.d_b(2))^2 - 1.0) - ...
                2.0*pi/obj.z(2))*(sqrt((obj.d_a(1)/obj.d_b(1))^2 - 1.0) - ...
                (obj.eps_alpha - 1.0)*2.0*pi/obj.z(1)));
            
            Z_BD = ones(1, 2)*nan;
            if((obj.eps_beta == 0.0) && (obj.eps_alpha > 1.0))
                if(M_1 > 1.0)
                    Z_BD(1) = M_1;
                else
                    Z_BD(1) = 1.0;
                end
                
                if(M_2 > 1.0)
                    Z_BD(2) = M_2;
                else
                    Z_BD(2) = 1.0;
                end
            elseif((obj.eps_alpha > 1.0) && (obj.eps_beta >= 1.0))
                Z_BD(1) = 1.0;
                Z_BD(2) = 1.0;
            elseif((obj.eps_alpha > 1.0) && (obj.eps_beta <  1.0))
                Z_BD(1) = M_1 - obj.eps_beta*(M_1 - 1.0);
                Z_BD(2) = M_2 - obj.eps_beta*(M_2 - 1.0);
            end
        end
        
        function [Z_L, Z_v] = lubrication_velocity_factor(~, sigma_Hlim, nu_40, v)
            if(sigma_Hlim  < 850.0) % [N/mm^2]
                C_ZL = 0.83;
            elseif((850.0 <= sigma_Hlim) && (sigma_Hlim  < 1200.0))
                C_ZL = sigma_Hlim/4375.0 + 0.6357;
            else
                C_ZL = 0.91;
            end
            
            % Lubricant factor:
            Z_L = C_ZL + 4.0*(1.0 - C_ZL)/(1.2 + 134.0/nu_40)^2;
            
            % Velocity factor:
            C_Zv = C_ZL + 0.02;
            Z_v = C_Zv + 2.0*(1.0 - C_Zv)/sqrt(0.8 + 32.0/v);
        end
        
        function Z_R = roughness_factor(obj, sigma_Hlim)
            rho_1 = 0.5*obj.d_b(1)*tand(obj.alpha_wt);
            rho_2 = 0.5*obj.d_b(2)*tand(obj.alpha_wt);
            
            rho_red = (rho_1*rho_2)/(rho_1 + rho_2);
            
            R_z10 = obj.R_z*power(10.0/rho_red, 1.0/3.0);
            
            if(sigma_Hlim < 850.0) % [N/mm^2]
                C_ZR = 0.15;
            elseif((850.0 <= sigma_Hlim) && (sigma_Hlim < 1200.0))
                C_ZR = 0.32 - sigma_Hlim*2.0e-4;
            else
                C_ZR = 0.08;
            end
            
            Z_R = power(3.0/R_z10, C_ZR);
        end
        
        function Z_NT = life_factor(~, line, N)
            switch line
                case 1
                    % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF (when limited pitting is permitted)
                    x = [6.0e5 1.0e7 1.0e9 1.0e10];
                    y = [1.6   1.3   1.0   0.85];
                case 2
                    % St, V, GGG (perl. bai.), GTS (perl.), Eh, IF
                    x = [1.0e5 5.0e7 1.0e9 1.0e10];
                    y = [1.6   1.0   1.0   0.85];
                case 3
                    % GG, GGG (ferr.), NT (nitr.), NV (nitr.)
                    x = [1.0e5 2.0e6 1.0e10];
                    y = [1.3   1.0   0.85];
                case 4
                    % NV (nitrocar.)
                    x = [1.0e5 2.0e6 1.0e10];
                    y = [1.1   1.0   0.85];
                otherwise
                    error("prog:input", "Invalid input [%d].\nValid options are 1 to 4.\n", line);
            end

            if(N <= x(1))
                Z_NT = y(1);
            elseif(N > x(end))
                Z_NT = y(end);
            else
                Z_NT = interp1(log(x), y, log(N));
            end
        end
    end
    
    %% Set methods:
    
    %% Get methods:
    methods
        function val = get.T_1(obj)
            val = (obj.P_rated*1.0e3)/(obj.n_out*pi/30.0);
        end
        
        function val = get.K_gamma(obj)
            Np = obj.N_p;
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
