classdef (Abstract) ISO_6336 < Gear_Set
    %ISO_6336 do some calculations following the ISO 6336 standard.
    % see also: KSISO_6336, MATISO_6336, Gear_Set
    %
    % References:
    %   [1] ISO 6336-1:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 1: Basic principles, introduction and general 
    % influence factors.
    %   [2] ISO 6336-2:2006 Calculation of load capacity of spur and 
    % helical gears -- Part 2: Calculation of surface durability (pitting).
    %   [3] ISO 6336-6:2019 Calculation of load capacity of spur and 
    % helical gears -- Part 6: Calculation of service life under variable 
    % load
    %   [4] ISO/TR 6336-30:2017 Calculation of load capacity of spur and 
    % helical gears -- Calculation examples for the application of ISO 6336 
    % parts 1, 2, 3, 5
    %   [5] Arnaudov, K., Karaivanov, D. (2019). Planetary Gear Trains. 
    % Boca Raton: CRC Press, https://doi.org/10.1201/9780429458521
    %
    
    properties(Access = public)
        P_rated;     % [kW],     Rated power
        S_Hmin;      % [-],      Minimum required safety factor for surface durability (pitting)
        S_Fmin;      % [-],      Minimum required safety factor for tooth bending strength
        L_h;         % [h],      Required life
        K_A;         % [-],      Application factor
        n_nominal;   % [1/min.], Nominal speed of each component(gears and carrier)
    end
    
    properties(Access = private)
        idx_fixed; % [-],        Index of the fixed element (if planetary)
        idx_input; % [-],        Index of the input element
        gear_set;
    end
    
    properties(Dependent, Access = public)
        speed_ratio;   % [-],      Speed ratio w.r.t. input element
        torque_ratio;  % [-],      Torque ratio w.r.t. input element
        T_nominal;     % [N-m],    Nominal torque applied on each element
        T_input;       % [N-m],    Input torque of the Gear_Set
        T_output;      % [N-m],    Output torque of the Gear_Set
        n_input;       % [1/min.], Input speed of the Gear_Set
        n_output;      % [1/min.], Output speed of the Gear_Set
        sigma_HP_ref;  % [N/mm^2], Permissive contact stress (reference)
        sigma_HP_stat; % [N/mm^2], Permissive contact stress (static)
        K_gamma;       % [-],      Mesh load factor
        F_tH;          % [N],      Determinant tangential load in a transverse plane for K_Halpha and K_Falpha
    end
    
    properties(Dependent, Access = private)
        idx_output; % [-],        Index of the output element
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
            
            default = {'P_rated'    , P_r, ...
                       'S_Hmin'     , 1.0, ...
                       'S_Fmin'     , 1.0, ...
                       'L_h'        , 50.0e3, ...
                       'K_A'        , 1.0, ...
                       'n_nominal'  , [n_1, NaN]};
            
            default = scaling_factor.process_varargin(default, varargin);
            
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
                         'R_a'          , gset.R_a, ...
                         'material'     , gset.material);
            
            obj.gear_set = gset;
            
            obj.P_rated = default.P_rated;
            obj.S_Hmin  = default.S_Hmin;
            obj.S_Fmin  = default.S_Fmin;
            obj.L_h     = default.L_h;
            obj.K_A     = default.K_A;

            [obj.n_nominal, ...
                obj.idx_fixed, obj.idx_input] = obj.set_gear_speed(default.n_nominal);
        end
        
        function [N, sigma] = Pitting_SN_curve(obj, varargin)
            N = logspace(1, 10, 100);
            default = {'N', N};
            default = scaling_factor.process_varargin(default, varargin);
            N = default.N;
            
            ZN = obj.life_factor(N);
            
            sigma = diag(obj.sigma_HP_ref)*ZN;
        end
        
        function disp(obj)
            % to be done...
        end
        
    end
    
    %% Calculation methods:
    methods
        function N_star = num_cycles_failure(obj)
            N_star = zeros(1, 2);
            sig_HPR = obj.sigma_HP_ref;
            
            for idx = 1:2
                func = @(x)(obj.sigma_H(idx).*obj.S_H(idx) - obj.life_factor(x).*sig_HPR(idx));
                x0 = obj.N_L(idx)*2.0;
                [N_star(idx),~] = fsolve(func, x0);
            end
        end
        
        function [num_cycle, sigmaH, torque_edges, speed_Max] = pitting_stress_bins(obj, torque, speed, time_step)
            %PITTING_STRESS_BINS performs LDD on the torque signal,
            % calculating the values of the pitting stress sigmaH for each
            % torque bin. The num_cycles is expanded to the required life
            % span L_h of the Gear_Set.
            
            [num_cycle, torque_edges, speed_Max] = ISO_6336.LDD(torque, speed, time_step);
            
            Max_power = torque_edges(2:end).*1.0e-3 .* ... % [kN-m]
                        speed_Max.*(pi/30.0);     % [rad/s] = [kW]
            
            obj_tmp = obj;
            obj_tmp.K_A = 1.0;
            
            sigmaH = zeros(2, length(speed_Max));
            
            for idx = 1:length(speed_Max)
                obj_tmp.P_rated = Max_power(idx);
                obj_tmp.n_nominal = obj_tmp.speed_ratio*abs(speed_Max(idx));
                sigmaH(:, idx) = obj_tmp.sigma_H';
            end
            
            time_duration = time_step*length(torque);
            num_cycle = num_cycle.*obj.L_h./time_duration;
            
            %% Example of plotting:
%             figure;
%             subplot(211)
%             histogram('binEdges' , num_cycle, ...
%                       'binCounts', sigmaH(1, :));
%             ylabel('\sigma_{H1}, [N/mm^2]');
%             
%             subplot(212)
%             histogram('binEdges' , num_cycle, ...
%                       'binCounts', sigmaH(2, :));
%             xlabel('N, [-]');
%             ylabel('\sigma_{H2}, [N/mm^2]');
%             
%             fig_axes = findobj(gcf, 'Type', 'Axes');
%             set(fig_axes, 'box', 'on');
%             set(fig_axes, 'xscale', 'log');
% 
        end
        
        %%
        function D = Simpack_damage_analysis(obj, data)
            field_name = fieldnames(data.load);
            
            data_tmp = data;
            
            if(strcmp(obj.configuration, 'planetary'))
                % Planet bearings:
                D_PL_A = zeros(1, obj.N_p);
                D_PL_B = D_PL_A;
                idx_PL  = find(contains(field_name, 'S_planet'));
                for idx = 1:obj.N_p
                    jdx = idx_PL(2*idx - 1);
                    data_tmp.load = data.load.(field_name{jdx});
                    D_PL_A(idx) = obj.bearing(1).Simpack_damage_analysis(data_tmp);
                    
                    jdx = idx_PL(2*idx);
                    data_tmp.load = data.load.(field_name{jdx});
                    D_PL_B(idx) = obj.bearing(2).Simpack_damage_analysis(data_tmp);
                end
                
                D_PL = [D_PL_A; D_PL_B];
                D_PL = reshape(D_PL, 1, 2*obj.N_p);
                
                % Planet carrier bearings:
                idx_PLC = find(contains(field_name, 'S_carrier'));
                jdx = idx_PLC(1);
                data_tmp.load = data.load.(field_name{jdx});
                D_PLC(1) = obj.bearing(3).Simpack_damage_analysis(data_tmp);
                
                jdx = idx_PLC(2);
                data_tmp.load = data.load.(field_name{jdx});
                D_PLC(2) = obj.bearing(4).Simpack_damage_analysis(data_tmp);
                
                D = [D_PL, D_PLC];
                
%                 idx_SP = find(contains(field_name, 'sun_planet'));
%                 name_field = fieldnames(data.speed);
%                 D_SP = zeros(1, obj.N_p);
%                 for idx = 1:obj.N_p
%                     jdx = idx_SP(idx);
%                     data_tmp.load = data.load.(field_name{jdx});
%                     
%                     [n, sigma_H, speed] = ...
%                     ISO_6336.LDD(data_tmp.load.ov_020.values, ...
%                                  data_tmp.speed.(name_field{1}).x.values, ...
%                                  data_tmp.time_step);
%                 end
                
%                 for idx = idx_SP
%                     ISO_6336.LDD(data.load.(field_name{idx}), ...
%                                  data.speed.(field_name{idx}), ...
%                                  data.time_step);
%                 end
%                 idx_RP = find(contains(field_name, '_ring_planet_'));
                
            else
                % Pinion bearings:
                idx_PIN  = find(contains(field_name, 'pinion_b'));
                D_pin = zeros(1, length(idx_PIN));
                for idx = 1:length(D_pin)
                    jdx = idx_PIN(idx);
                    data_tmp.load = data.load.(field_name{jdx});
                    D_pin(idx) = obj.bearing(idx).Simpack_damage_analysis(data_tmp);
                end
                
                % Wheel bearings:
                idx_WHE  = find(contains(field_name, 'wheel_b'));
                D_whe = zeros(1, length(idx_WHE));
                for idx = 1:length(D_whe)
                    jdx = idx_WHE(idx);
                    data_tmp.load = data.load.(field_name{jdx});
                    kdx = length(D_pin) + idx;
                    D_whe(idx) = obj.bearing(kdx).Simpack_damage_analysis(data_tmp);
                end
                
                D = [D_pin, D_whe];
                
%                 idx_PW = find(contains(field_name, 'pinion_wheel'));
%                 ISO_6336.LDD(data.load.(field_name{idx_PW}), ...
%                              data.speed.(field_name{idx_PW}), ...
%                              data.time_step);
            end
            
        end
        
        function show_LDD(obj, data)
            field_name = fieldnames(data.load);
            name_field = fieldnames(data.speed);
            
            data_idx = data;
            
            if(strcmp(obj.configuration, 'planetary'))
                % Planet bearings:
                figure;
                idx_PL  = find(contains(field_name, 'S_planet'));
                for idx = 1:obj.N_p
                    jdx = idx_PL(2*idx - 1);
                    data_idx.load = data.load.(field_name{jdx});
                    subplot(3, 2, 1);
                    hold on;
                    obj.bearing(1).show_LDD(data_idx);
                    
                    jdx = idx_PL(2*idx);
                    data_idx.load = data.load.(field_name{jdx});
                    subplot(3, 2, 3);
                    hold on;
                    obj.bearing(2).show_LDD(data_idx);
                end
                
                subplot(3, 2, 1);
                title('Planet bearing A');
                legend({'A', 'B', 'C'}, 'location', 'best');
                
                subplot(3, 2, 3);
                title('Planet bearing B');
                
                
                % Planet carrier bearings:
                idx_PLC = find(contains(field_name, 'S_carrier'));
                for idx = 1:length(idx_PLC)
                    jdx = idx_PLC(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    subplot(3, 2, 5);
                    hold on;
                    obj.bearing(idx + 2).show_LDD(data_idx);
                end
                title('Planet carrier bearing');
                
                % Sun-planet gear contact:
                jdx_SP = find(contains(field_name, 'sun_planet'));
                kdx_SP = find(contains(name_field, 'planet'));
                subplot(3, 2, 2)
                hold on;
                for idx = 1:obj.N_p
                    jdx = jdx_SP(idx);
                    kdx = kdx_SP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    [N, sigma_H] = ...
                    ISO_6336.LDD(data_idx.load.ov_020.values, ...
                                 data_idx.speed.(name_field{kdx}).x.values, ...
                                 data_idx.time_step);
                    
                    histogram('binEdges' , N, ...
                              'binCounts', sigma_H(2:end)*1.0e-6);
                    box on;
                end
                title('Max. Sun-Planet contact stress');
                
                % Ring-planet gear contact:
                jdx_RP = find(contains(field_name, 'ring_planet'));
                kdx_RP = find(contains(name_field, 'planet'));
                subplot(3, 2, 4)
                hold on;
                for idx = 1:obj.N_p
                    jdx = jdx_RP(idx);
                    kdx = kdx_RP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    [N, sigma_H] = ...
                    ISO_6336.LDD(data_idx.load.ov_020.values, ...
                                 data_idx.speed.(name_field{kdx}).x.values, ...
                                 data_idx.time_step);
                    
                    histogram('binEdges' , N, ...
                              'binCounts', sigma_H(2:end)*1.0e-6);
                end
                title('Max. Ring-Planet contact stress');
                
            else
                % Pinion bearings:
                figure;
                idx_PIN  = find(contains(field_name, 'pinion_b'));
                n = length(idx_PIN);
                subplot(3, 1, 1)
                hold on;
                for idx = 1:n
                    jdx = idx_PIN(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    obj.bearing(idx).show_LDD(data_idx);
                end
                title('Pinion bearing');
                
                % Wheel bearings:
                idx_WHE  = find(contains(field_name, 'wheel_b'));
                n = length(idx_WHE);
                subplot(3, 1, 2)
                hold on;
                for idx = 1:n
                    jdx = idx_WHE(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    kdx = n + idx;
                    
                    obj.bearing(kdx).show_LDD(data_idx);
                end
                title('Wheel bearing');
                
                % Pinion-Wheel gear contact:
                jdx_PW = contains(field_name, 'pinion_wheel');
                kdx_PW = contains(name_field, 'pinion');
                name_field = fieldnames(data.speed);
                
                [N, sigma_H] = ISO_6336.LDD(data.load. (field_name{jdx_PW}).ov_020.values, ...
                                            data.speed.(name_field{kdx_PW}).x.values, ...
                                            data.time_step);

%                 figure;
                subplot(3, 1, 3)
                hold on;
                histogram('binEdges' , N, ...
                          'binCounts', sigma_H(2:end)*1.0e-6);
                title('Max. Pinion-Wheel contact stress');
                                        
            end
            
            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'box', 'on');
            set(fig_axes, 'xScale', 'log');

        end
        
    end
    
    methods(Static)
        function [num_cycle, load_edges, speed_Max, bin] = LDD(load_signal, load_speed, time_step)
            %LDD produces the load duration distribution of a load_signal
            % with rotational speed load_speed in [1/min.]. time_step is 
            % the inverse of the sampling frequency of both signals.
            %
            % References:
            % [1] A. R. Nejad, Z. Gao, and T. Moan, "On long-term fatigue 
            % damage and reliability analysis of gears under wind loads in 
            % offshore wind turbine drivetrains" Int. J. Fatigue, vol. 61,
            % pp. 116-128, Apr. 2014. 10.1016/j.ijfatigue.2013.11.023
            %

            load_signal = abs(load_signal);
            load_speed = abs(load_speed);
            
            [load_min, load_Max] = bounds(load_signal);
            num_bins = 101; % according to [1], 4.6
            bin_size = (load_Max - load_min)/num_bins;
            if(load_min > bin_size)
                load_min = load_min - bin_size;
            end
            load_Max = load_Max + bin_size;
            load_edges = linspace(load_min, load_Max, num_bins);
            
            % [N, edges, bin] = histcounts(signal, nbins OR edges)
            [N, load_edges, bin] = histcounts(load_signal, load_edges);
            
            speed_Max = zeros(1, num_bins - 1);
            for idx = 1:(num_bins - 1)
                if(N(idx) ~= 0.0)
                    range = (load_signal >= load_edges(idx)) & ...
                            (load_signal <= load_edges(idx + 1));
                    speed_Max(idx) = max(load_speed(range));
                end
            end
            
            % descending order:
            load_edges = fliplr(load_edges);
            speed_Max  = fliplr(speed_Max);
            N          = fliplr(N);
            
            % Assuming that each load is constant during a time step:
            load_time = N*time_step;
            
            % according to [1], (10):
            num_cycle = load_time.*speed_Max/60.0;
            num_cycle = [0.0, cumsum(num_cycle)];
            
            %% Example of plotting:
%             figure;
%             subplot(211)
%             histogram('binEdges', fliplr(load_edges), ...
%                       'binCounts', fliplr(N));
%             title('Histogram');
%             xlabel('Load, [NA]');
%             ylabel('Num. cycles, [-]');
%             
%             subplot(212)
%             histogram('binEdges' , num_cycle, ...
%                       'binCounts', load_edges(2:end));    ...
%             title('Load duration distribution - LDD');
%             xlabel('Cum. num. cycles, [-]');
%             ylabel('Load, [NA]');
%             set(gca, 'xscale', 'log')
%             
%             fig_axes = findobj(gcf, 'Type', 'Axes');
%             set(fig_axes, 'box', 'on');
%             set(fig_axes, 'xscale', 'log');
% 
        end
        
    end

    methods(Access = private)
        function val = F_t_spectra(obj, torque)
            val = 2.0e3*(torque/obj.d(1))/obj.N_p;
            val = abs(val);
        end
        
        function val = F_tH_spectra(obj, torque)
            val = 2.0e3*(torque/obj.d(1))/obj.N_p;
            val = abs(val);
        end
        
        function val = dynamic_factor(obj)
            z1 = obj.z(1);
            uu = obj.u;
            cond = (obj.v_pitch_line*z1/100.0)*sqrt((uu.^2)/(1.0 + uu.^2));
            if(cond < 3.0) % [m/s]
                warning('ISO_6336:KV', ['Calculating K_v using method B', ...
                    ' outside of its useful range. ', ...
                'More info at the end of Sec. 6.3.2 of ISO 6336-1.']);
            end
            
            rho = [obj.material.rho];
            
            % transverse base pitch deviation, [um]:
            fpb = obj.f_pb;
            
            % profile form deviation, [um]
            f_falpha = max(obj.f_falpha);
            
            % Estimated running allowance (pitch deviation):
            y_p = obj.y_alpha;
            
            % Estimated running allowance (flank deviation):
            y_f = f_falpha*75.0e-3;
            
            % transverse effective base pitch deviation after running-in, [um]:
            f_pbeff = fpb - y_p;
            
            % effective profile form deviation after running-in, [um]
            f_falphaeff = f_falpha - y_f; % 0.925*y_f
            
            if(strcmp(obj.configuration, 'parallel'))
                % Based on Sec. 6.5.9 of ISO 6336-1, Eq. (30), assuming gears:
                % - of solid construction, and
                % - with the same density
                
                num = pi*rho(1)*(obj.u^2)*obj.d_m(1)^4;
                den = 8.0*(obj.u^2 + 1.0)*obj.d_b(1)^2;
                m_red = num/den;
                
                % Resonance running speed, Eq. (6) [1]:
                n_E1 = 3.0e4*sqrt(obj.c_gamma_alpha/m_red)/(pi*obj.z(1)); % [1/min]
                
                % Resonance ratio, Eq. (9) [1]:
                N = obj.n_output/n_E1;
                
                val = obj.dynamic_factor_from_range(N, obj.C_a, ...
                    f_pbeff, f_falphaeff);
            elseif(strcmp(obj.configuration, 'planetary'))
                m_sun = rho(1)*(pi/8.0)*(obj.d_m(1)^4)/(obj.d_b(1)^2); % sun
                m_pla = rho(2)*(pi/8.0)*(obj.d_m(2)^4)/(obj.d_b(2)^2); % planet
                m_red1 = (m_pla*m_sun)/(m_sun + m_pla*obj.N_p);
                m_red2 =  m_pla;
                
                % Resonance running speed:
                n_E11 = 3.0e4*sqrt(obj.cprime/m_red1)*1.0/(pi*obj.z(1)); % [1/min]
                n_E12 = 3.0e4*sqrt(obj.cprime/m_red2)*1.0/(pi*obj.z(2)); % [1/min]
                
                % Resonance ratio:
                N_1 =  obj.n_output/n_E11;
                N_2 = (obj.n_output/obj.u)/n_E12;

                K_v1 = obj.dynamic_factor_from_range(N_1, obj.C_a, ...
                    f_pbeff, f_falphaeff);
                
                K_v2 = obj.dynamic_factor_from_range(N_2, obj.C_a, ...
                    f_pbeff, f_falphaeff);
                
                val = max(K_v1, K_v2);
            end
        end
        
        function K_v = dynamic_factor_from_range(obj, N, C_a, f_pbeff, f_falphaeff)
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
            B_p = cp*f_pbeff/obj.line_load;        % Eq. (15)
            B_f = cp*f_falphaeff/obj.line_load;    % Eq. (16)
            B_k = abs(1.0 - cp*C_a/obj.line_load); % Eq. (17)
            
            % Dynamic factor:
            if(obj.line_load < 100) % [N/mm]
                N_S = 0.5 + 0.35*sqrt(obj.line_load/100.0); % Eq. (11), [1]
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
        
        function [Z_L, Z_v] = lubrication_velocity_factor(obj)

            % [N/mm^2],  Allowable contact stress number:
            sig_Hlim = min(obj.sigma_Hlim);
            
            if(sig_Hlim  < 850.0) % [N/mm^2]
                C_ZL = 0.83;
            elseif((850.0 <= sig_Hlim) && (sig_Hlim  < 1200.0))
                C_ZL = sig_Hlim/4375.0 + 0.6357;
            else
                C_ZL = 0.91;
            end
            
            % Lubricant factor:
            Z_L = C_ZL + 4.0*(1.0 - C_ZL)/(1.2 + 134.0/obj.nu_40)^2;
            
            % Velocity factor:
            C_Zv = C_ZL + 0.02;
            Z_v = C_Zv + 2.0*(1.0 - C_Zv)/sqrt(0.8 + 32.0/obj.v_pitch_line);
        end
        
    end
    
    %% Set methods:
    methods
    end
    
    %% Get methods:
    methods
        function val = get.speed_ratio(obj)
            val = obj.n_nominal./obj.n_input;
        end
        
        function val = get.torque_ratio(obj)
            val = obj.T_nominal./obj.T_input;
        end
        
        function val = get.F_tH(obj)
            val = obj.F_t*obj.K_gamma* ...
                          obj.K_A    * ...
                          obj.K_v    * ...
                          obj.K_Hbeta;
        end
        
        function val = get.sigma_HP_ref(obj)
            ZNT = 1.0;
            val = [obj.material(1:2).sigma_Hlim].*ZNT*obj.Z_L* ...
                                                      obj.Z_v* ...
                                                      obj.Z_R* ...
                                                      obj.Z_W* ...
                                                      obj.Z_X./obj.S_Hmin;
        end
        
        function val = get.sigma_HP_stat(obj)
            val = [obj.material(1:2).sigma_Hlim].*obj.Z_NT* ...
                                                  obj.Z_L* ...
                                                  obj.Z_v* ...
                                                  obj.Z_R* ...
                                                  obj.Z_W* ...
                                                  obj.Z_X./obj.S_Hmin;
        end
        
        function val = get.T_nominal(obj)
            val = (obj.P_rated*1.0e3)./(obj.n_nominal*pi/30.0);
            val(2) = val(2)/obj.N_p;
            val(val == inf) = nan;
        end
        
        function val = get.T_input(obj)
            val = obj.T_nominal(obj.idx_input);
        end
        
        function val = get.T_output(obj)
            val = obj.T_nominal(obj.idx_output);
        end
        
        function val = get.n_input(obj)
            val = obj.n_nominal(obj.idx_input);
        end
        
        function val = get.n_output(obj)
            val = obj.n_nominal(obj.idx_output);
        end
        
        function val = get.idx_output(obj)
            if(strcmpi(obj.configuration, 'parallel'))
                val = setdiff([1 2], obj.idx_input);
            elseif(strcmpi(obj.configuration, 'planetary'))
                switch(obj.idx_fixed)
                    case 1
                        val = setdiff([3 4], obj.idx_input);
                    case 2
                        error('ISO_6336:idx_out', 'Direct drive not implemented yet.');
                    case 3
                        val = setdiff([1 4], obj.idx_input);
                    case 4
                        val = setdiff([1 3], obj.idx_input);
                    otherwise
                        error('ISO_6336:idx_out', 'Invalid index.');
                end
            end
        end
        
        function val = get.K_gamma(obj)
            switch(obj.N_p)
                case 3
                    val = 1.1;
                case 4
                    val = 1.25;
                case 5
                    val = 1.35;
                case 6
                    val = 1.44;
                case 7
                    val = 1.47;
                otherwise
                    val = 1.0;
            end
        end
        
    end
end
