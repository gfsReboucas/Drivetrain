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
%         gear_set;
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
        J_eq;          % [kg-m^2], Equivalent mass moment of inertia
        k_T;           % [-],      Scaling factor for the input torque 
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
        function obj_sca = scale_by(obj_ref, gamma)
            if(~isa(gamma, 'scaling_factor'))
                error('gamma should be a [SCALING_FACTOR] object.');
            end
                                                        % Scaling factors for:
            gamma_one = scaling_factor({'m_n', 1.0, ... % normal module,         [mm]
                                        'b'  , 1.0, ... % face width,            [mm]
                                        'd'  , 1.0, ... % output shaft diameter, [mm]
                                        'L'  , 1.0});   % output shaft length,   [mm]

            gamma = gamma_one.update(gamma);
            
            m_n_sca = Rack.module(obj_ref.m_n*gamma('m_n'), 'round_0125', '');
            
            gamma('m_n') = m_n_sca/obj_ref.m_n;
            
            ref_shaft = obj_ref.output_shaft;
            shaft_sca = Shaft('d'       , ref_shaft.d*gamma('d'), ...
                              'L'       , ref_shaft.L*gamma('L'), ...
                              'bearing' , ref_shaft.bearing, ...
                              'material', ref_shaft.material);
                          
            gset = Gear_Set('configuration', obj_ref.configuration, ...
                            'm_n'          , obj_ref.m_n*gamma('m_n'), ...
                            'alpha_n'      , obj_ref.alpha_n, ...
                            'z'            , obj_ref.z, ...
                            'b'            , obj_ref.b*gamma('b'), ...
                            'x'            , obj_ref.x, ...
                            'beta'         , obj_ref.beta, ...
                            'k'            , obj_ref.k, ...
                            'bore_ratio'   , obj_ref.bore_ratio, ...
                            'N_p'          , obj_ref.N_p, ...
                            'a_w'          , obj_ref.a_w*gamma('m_n'), ...
                            'rack_type'    , obj_ref.type, ...
                            'bearing'      , obj_ref.bearing, ...
                            'shaft'        , shaft_sca, ...
                            'Q'            , obj_ref.Q, ...
                            'R_a'          , obj_ref.R_a, ...
                            'material'     , obj_ref.material);
            
            if(strcmpi(obj_ref.configuration, 'parallel'))
                n_nom = [1.0, nan];
            elseif(strcmpi(obj_ref.configuration, 'planetary'))
                n_nom = [nan, nan, 0.0, 1.0];
            end
            n_nom = n_nom*obj_ref.n_input;
            
            constructor = str2func(class(obj_ref));
            obj_sca = constructor(gset, 'P_rated'  , obj_ref.P_rated, ...
                                        'S_Hmin'   , obj_ref.S_Hmin, ...
                                        'S_Fmin'   , obj_ref.S_Fmin, ...
                                        'L_h'      , obj_ref.L_h, ...
                                        'K_A'      , obj_ref.K_A, ...
                                        'n_nominal', n_nom);
        end
        
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
%                       'binCounts', sigmaH(1, :), ...
%                       'displayStyle', 'stairs', ...
%                       'lineWidth'   , 2.0);
%             ylabel('\sigma_{H1}, [N/mm^2]');
%             
%             subplot(212)
%             histogram('binEdges' , num_cycle, ...
%                       'binCounts', sigmaH(2, :), ...
%                       'displayStyle', 'stairs', ...
%                       'lineWidth'   , 3);
%             xlabel('N, [-]');
%             ylabel('\sigma_{H2}, [N/mm^2]');
%             
%             fig_axes = findobj(gcf, 'Type', 'Axes');
%             set(fig_axes, 'box', 'on');
%             set(fig_axes, 'xscale', 'log');
% 
        end
        
        %%
        function load_spectrum_KISSsoft(obj, file_name, torque_signal, load_speed, time_step)
            [~, torque_edges, speed_edges, N] = ISO_6336.LDD(torque_signal, load_speed, time_step);
            
            power_edges = torque_edges(2:end).*1.0e-3 .* ... % [kN-m]
                           speed_edges.*(pi/30.0);    % [rad/s] = [kW]

            freq_factor  = N./sum(N);
            power_factor = power_edges./obj.P_rated;
            speed_factor = speed_edges./obj.n_input;

            tab = table(freq_factor',  power_factor',  speed_factor', 'variableNames', ...
                     {'#frequency'  , 'power'       , 'speed'});
            
            writetable(tab, sprintf('%s.dat', file_name));
        end
        
        function [N_L, sigH, sigF, par] = read_SNCurve_KS(obj, file_name)
            %READ_SNCURVE_KS reads the S-N curve obtained using KISSsoft
            % software. To generate this curve in KISSsoft go to the menu
            % Graphics > Evaluation > Woehler lines (S-N curves). The data
            % from these curves can be found at %temp%\KISS_XY\Z01.H1.TMP,
            % where %temp% is Windows temporary folder.
            %
            
            arguments
                obj;
                file_name = [];
            end
            
            if(strcmpi(obj.configuration, 'parallel'))
                format_spec = '  -        %f     %f N/mm²    %f N/mm²      %f N/mm² (2)  (2)     %f N/mm² (3)  (3)  ';
                idx_F = 2:3;
                idx_H = idx_F + 2;
            elseif(strcmpi(obj.configuration, 'planetary'))
                format_spec = '  -        %f     %f N/mm²    %f N/mm²    %f N/mm²      %f N/mm² (3)  (3)     %f N/mm² (4)  (4)     %f N/mm² (5)  (5)  ';
                idx_F = 2:4;
                idx_H = idx_F + 3;
            end
            
            if(isempty(file_name))
                file_name = tempdir;
                file_name = strrep(file_name, ':', '');
                idx = find(ismember(file_name, 'Temp'), 1, 'last') + 1;
                file_name = ['\\tsclient\', file_name(1:idx), 'KISS_*\Z01-H1.TMP'];
            end
            
            file = dir(file_name);
            file_full = [file.folder, '\', file.name];
            file_ID = fopen(file_full, 'r');
            
            line_start = 7;
            line_start = line_start - 1;
            
            C = textscan(file_ID, format_spec, 1000, 'delimiter', '\n', 'headerlines', line_start);
            N_L = C{1};
            sigH = cell2mat(C(idx_H));
            sigF = cell2mat(C(idx_F));
            
            par = struct;
            for idx = 1:numel(obj.z)
                %% Root-bending:
                jdx = find(sigF(:, idx) <  sigF(1  , idx), 1, 'first');
                kdx = find(sigF(:, idx) <= sigF(end, idx), 1, 'first');
                
                b = nlinfit(log([N_L(jdx:kdx), sigF(jdx:kdx, idx)]), ... % x
                            zeros(kdx - jdx + 1, 1), ... % y
                            @log_SNC, ones(2, 1));
                par.F.m(idx) = b(1);
                par.F.logK(idx) = b(2);

                %% Pitting:
                jdx = find(sigH(:, idx) <  sigH(1  , idx), 1, 'first');
                kdx = find(sigH(:, idx) <= sigH(end, idx), 1, 'first');
                
                b = nlinfit(log([N_L(jdx:kdx), sigH(jdx:kdx, idx)]), ... % x
                            zeros(kdx - jdx + 1, 1), ... % y
                            @log_SNC, ones(2, 1));
                par.H.m(idx) = b(1);
                par.H.logK(idx) = b(2);
            end
            
            function e = log_SNC(b, x)
                % Originally K = N S ^ m, then, after isolating S and
                % applying log on both sides one obtains: b_1 y = b_2 - x,
                % where: x = log(N), y = log(S), b_1 = m, b_2 = log(K).
                % 
                e = b(1)*x(:, 2) + x(:, 1) - b(2);
            end
        end
        
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
        
        function show_LDD(obj, data, name)
            field_name = fieldnames(data.load);
            name_field = fieldnames(data.speed);
            
            data_idx = data;
            
            if(strcmp(obj.configuration, 'planetary'))
                % Planet bearings:
                figure('units', 'centimeters', 'position', [5.0, 5.0, 75.0 15.0]);
                idx_PL  = find(contains(field_name, 'S_planet'));
                for idx = 1:obj.N_p
                    jdx = idx_PL(2*idx - 1);
                    data_idx.load = data.load.(field_name{jdx});
                    subplot(3, 5, 5*idx - 4);
                    obj.bearing(1).show_LDD(data_idx);
                    ylabel('P, [N]')
                    title(sprintf('%s.%d', obj.bearing(1).name, idx));
                    
                    jdx = idx_PL(2*idx);
                    data_idx.load = data.load.(field_name{jdx});
                    subplot(3, 5, 5*idx - 3);
                    obj.bearing(2).show_LDD(data_idx);
                    ylabel('P, [N]')
                    title(sprintf('%s.%d', obj.bearing(2).name, idx));
                end
                subplot(3, 5, 11)
                xlabel('N, [-]')
                
                subplot(3, 5, 12)
                xlabel('N, [-]')
                
                % Planet carrier bearings:
                idx_PLC = find(contains(field_name, 'S_carrier'));
                for idx = 1:length(idx_PLC)
                    jdx = idx_PLC(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    subplot(3, 5, 5*idx - 2);
                    hold on;
                    obj.bearing(idx + 2).show_LDD(data_idx);
                    ylabel('P, [N]')
                    title(sprintf('%s.%d', obj.bearing(idx + 2).name, idx));
                end
                xlabel('N, [-]')
                
                % Sun-planet gear contact:
                jdx_SP = find(contains(field_name, 'sun_planet'));
                kdx_SP = find(contains(name_field, 'planet'));
                for idx = 1:obj.N_p
                    jdx = jdx_SP(idx);
                    kdx = kdx_SP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    [N, sigma_H] = ...
                    ISO_6336.LDD(data_idx.load.ov_020.values, ...
                                 data_idx.speed.(name_field{kdx}).x.values, ...
                                 data_idx.time_step);
                             
                    subplot(3, 5, 5*idx - 1)
                    histogram('binEdges'    , N, ...
                              'binCounts'   , sigma_H(2:end)*1.0e-6, ...
                              'displayStyle', 'stairs', ...
                              'lineWidth'   , 2.0);
                    ylabel('\sigma_H, [N/mm^2]')
                end
                xlabel('N, [-]')
                subplot(3, 5, 4)
                title('Max. Sun-Planet contact stress');
                
                % Ring-planet gear contact:
                jdx_RP = find(contains(field_name, 'ring_planet'));
                kdx_RP = find(contains(name_field, 'planet'));
                for idx = 1:obj.N_p
                    jdx = jdx_RP(idx);
                    kdx = kdx_RP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    [N, sigma_H] = ...
                    ISO_6336.LDD(data_idx.load.ov_020.values, ...
                                 data_idx.speed.(name_field{kdx}).x.values, ...
                                 data_idx.time_step);
                    
                    subplot(3, 5, 5*idx)
                    histogram('binEdges'    , N, ...
                              'binCounts'   , sigma_H(2:end)*1.0e-6, ...
                              'displayStyle', 'stairs', ...
                              'lineWidth'   , 2.0);
                    ylabel('\sigma_H, [N/mm^2]')
                end
                xlabel('N, [-]')
                subplot(3, 5, 5)
                title('Max. Ring-Planet contact stress');
                
            else
                % Pinion bearings:
                figure('units', 'centimeters', 'position', [5.0, 5.0, 45.0 15.0]);
                idx_PIN  = find(contains(field_name, 'pinion_b'));
                n_PIN = length(idx_PIN);
                for idx = 1:n_PIN
                    jdx = idx_PIN(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    subplot(3, 3, 3*idx - 2)
                    obj.bearing(idx).show_LDD(data_idx);
                    ylabel('P, [N]')
                    title(sprintf('%s.%d', obj.bearing(idx).name, idx));
                end
                
                % Wheel bearings:
                idx_WHE  = find(contains(field_name, 'wheel_b'));
                n_WHE = length(idx_WHE);
                for idx = 1:n_WHE
                    jdx = idx_WHE(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    kdx = n_PIN + idx;
                    
                    subplot(3, 3, 3*idx - 1)
                    obj.bearing(kdx).show_LDD(data_idx);
                    ylabel('P, [N]')
                    title(sprintf('%s.%d', obj.bearing(kdx).name, idx));
                end
                
                % Pinion-Wheel gear contact:
                jdx_PW = contains(field_name, 'pinion_wheel');
                kdx_PW = contains(name_field, 'pinion');
                name_field = fieldnames(data.speed);
                
                [N, sigma_H] = ISO_6336.LDD(data.load. (field_name{jdx_PW}).ov_020.values, ...
                                            data.speed.(name_field{kdx_PW}).x.values, ...
                                            data.time_step);

                subplot(3, 3, 3:3:9)
                histogram('binEdges'    , N, ...
                          'binCounts'   , sigma_H(2:end)*1.0e-6, ...
                          'displayStyle', 'stairs', ...
                          'lineWidth'   , 2.0);
                title('Max. Pinion-Wheel contact stress');
                                        
            end
            
            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'box', 'on');
            set(fig_axes, 'xScale', 'log');
            
            fig_axes(end).Title.String = [name, ' - ', fig_axes(end).Title.String];

        end
        
        function show_histogram(obj, data, name)
            field_name = fieldnames(data.load);
            name_field = fieldnames(data.speed);
            
            data_idx = data;
            
            figure('units', 'centimeters', 'position', [5.0, 5.0, 50.0 20.0]);
            if(strcmp(obj.configuration, 'planetary'))
                % Planet bearings:
                idx_PL  = find(contains(field_name, 'S_planet'));
                for idx = 1:obj.N_p
                    jdx = idx_PL(2*idx - 1);
                    data_idx.load = data.load.(field_name{jdx});
                    subplot(3, 5, 5*idx - 4);
                    obj.bearing(1).show_histogram(data_idx);
                    ylabel('N, [-]')
                    title(sprintf('%s.%d', obj.bearing(1).name, idx));
                    
                    jdx = idx_PL(2*idx);
                    data_idx.load = data.load.(field_name{jdx});
                    subplot(3, 5, 5*idx - 3);
                    obj.bearing(2).show_histogram(data_idx);
                    ylabel('N, [-]')
                    title(sprintf('%s.%d', obj.bearing(2).name, idx));
                end
                subplot(3, 5, 11)
                xlabel('P, [N]')
                
                subplot(3, 5, 12)
                xlabel('P, [N]')
                
                % Planet carrier bearings:
                idx_PLC = find(contains(field_name, 'S_carrier'));
                for idx = 1:length(idx_PLC)
                    jdx = idx_PLC(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    subplot(3, 5, 5*idx - 2);
                    hold on;
                    obj.bearing(idx + 2).show_histogram(data_idx);
                    ylabel('N, [-]')
                    title(sprintf('%s.%d', obj.bearing(idx + 2).name, idx));
                end
                xlabel('P, [N]')
                
                % Sun-planet gear contact:
                jdx_SP = find(contains(field_name, 'sun_planet'));
                kdx_SP = find(contains(name_field, 'planet'));
                for idx = 1:obj.N_p
                    jdx = jdx_SP(idx);
                    kdx = kdx_SP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    [N, sigma_H, speed] = ...
                    ISO_6336.LDD(data_idx.load.ov_020.values, ...
                                 data_idx.speed.(name_field{kdx}).x.values, ...
                                 data_idx.time_step);
                    N = diff(N)*60.0./(speed.*data.time_step);
                    N(isnan(N)) = 0;
            
                    subplot(3, 5, 5*idx - 1)
                    histogram('binEdges'    , fliplr(sigma_H)*1.0e-6, ...
                              'binCounts'   , fliplr(N), ...
                              'displayStyle', 'stairs', ...
                              'lineWidth'   , 2.0);
                    ylabel('N, [-]')
                end
                xlabel('\sigma_H, [N/mm^2]');

                subplot(3, 5, 4)
                title('Max. Sun-Planet contact stress');
                
                % Ring-planet gear contact:
                jdx_RP = find(contains(field_name, 'ring_planet'));
                kdx_RP = find(contains(name_field, 'planet'));
                for idx = 1:obj.N_p
                    jdx = jdx_RP(idx);
                    kdx = kdx_RP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    [N, sigma_H, speed] = ...
                    ISO_6336.LDD(data_idx.load.ov_020.values, ...
                                 data_idx.speed.(name_field{kdx}).x.values, ...
                                 data_idx.time_step);
                    N = diff(N)*60.0./(speed.*data.time_step);
                    N(isnan(N)) = 0;

                    subplot(3, 5, 5*idx)
                    histogram('binEdges'    , fliplr(sigma_H)*1.0e-6, ...
                              'binCounts'   , fliplr(N), ...
                              'displayStyle', 'stairs', ...
                              'lineWidth'   , 2.0);
                    ylabel('N, [-]')
                end
                xlabel('\sigma_H, [N/mm^2]');
                
                subplot(3, 5, 5)
                title('Max. Ring-Planet contact stress');
                
            else
                % Pinion bearings:
                idx_PIN  = find(contains(field_name, 'pinion_b'));
                n_PIN = length(idx_PIN);
                for idx = 1:n_PIN
                    jdx = idx_PIN(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    subplot(3, 3, 3*idx - 2)
                    obj.bearing(idx).show_histogram(data_idx);
                    ylabel('N, [-]')
                    title(sprintf('%s.%d', obj.bearing(idx).name, idx));
                end
                xlabel('P, [N]');
                
                % Wheel bearings:
                idx_WHE  = find(contains(field_name, 'wheel_b'));
                n_WHE = length(idx_WHE);
                for idx = 1:n_WHE
                    jdx = idx_WHE(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    kdx = n_PIN + idx;
                    
                    subplot(3, 3, 3*idx - 1)
                    obj.bearing(kdx).show_histogram(data_idx);
                    title(sprintf('%s.%d', obj.bearing(kdx).name, idx));
                end
                xlabel('P, [N]');
                
                % Pinion-Wheel gear contact:
                jdx_PW = contains(field_name, 'pinion_wheel');
                kdx_PW = contains(name_field, 'pinion');
                name_field = fieldnames(data.speed);
                
                [N, sigma_H, speed] = ISO_6336.LDD(data.load. (field_name{jdx_PW}).ov_020.values, ...
                                                   data.speed.(name_field{kdx_PW}).x.values, ...
                                                   data.time_step);
                N = diff(N)*60.0./(speed.*data.time_step);
                N(isnan(N)) = 0;

                subplot(3, 3, 3:3:9)
                histogram('binEdges'    , fliplr(sigma_H)*1.0e-6, ...
                          'binCounts'   , fliplr(N), ...
                          'displayStyle', 'stairs', ...
                          'lineWidth'   , 2.0);
                ylabel('N, [-]')
                title('Max. Pinion-Wheel contact stress');
            end
            
            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'box', 'on');
%             set(fig_axes, 'xScale', 'log');
            
            fig_axes(end).Title.String = [name, ' - ', fig_axes(end).Title.String];

        end
        
        function [a, b] = Weibull(obj, data, name)
            field_name = fieldnames(data.load);
            
            data_idx = data;
            
            if(strcmp(obj.configuration, 'planetary'))
                % Planet bearings:
                idx_PL  = find(contains(field_name, 'S_planet'));
                a_PL = zeros(size(idx_PL));
                b_PL = zeros(size(idx_PL));
                
                for idx = 1:obj.N_p
                    jdx = idx_PL(2*idx - 1);
                    data_idx.load = data.load.(field_name{jdx});

                    [a, b] = obj.bearing(1).Weibull(data_idx);
                    fprintf('%s\t%s%d\ta = %e\tb = %e\n', name, obj.bearing(1).name, ...
                                                        idx, a, b);

                    a_PL(2*idx - 1) = a;
                    b_PL(2*idx - 1) = b;
                    
                    jdx = idx_PL(2*idx);
                    data_idx.load = data.load.(field_name{jdx});

                    [a, b] = obj.bearing(2).Weibull(data_idx);
                    fprintf('%s\t%s%d\ta = %e\tb = %e\n', name, obj.bearing(2).name, ...
                                                        idx, a, b);

                    a_PL(2*idx) = a;
                    b_PL(2*idx) = b;
                end
                
                % Planet carrier bearings:
                idx_PLC = find(contains(field_name, 'S_carrier'));
                a_PLC = zeros(size(idx_PLC));
                b_PLC = zeros(size(idx_PLC));
                for idx = 1:length(idx_PLC)
                    jdx = idx_PLC(idx);
                    data_idx.load = data.load.(field_name{jdx});

                    [a, b] = obj.bearing(idx + 2).Weibull(data_idx);
                    fprintf('%s\t%s%d\ta = %e\tb = %e\n', name, obj.bearing(idx + 2).name, ...
                                                        idx, a, b);

                    a_PLC(idx) = a;
                    b_PLC(idx) = b;
                end
                
                % Sun-planet gear contact:
                jdx_SP = find(contains(field_name, 'sun_planet'));
                a_SP = zeros(obj.N_p, 1);
                b_SP = zeros(obj.N_p, 1);
                for idx = 1:obj.N_p
                    jdx = jdx_SP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                                 
                    wbull = fitdist(data_idx.load.ov_020.values + eps, 'weibull');
                    
                    a = wbull.a; % scale param.
                    b = wbull.b; % shape param.
                    
                    fprintf('%s:\tSP%d\ta = %e\tb = %e\n', name, idx, a, b);
                    a_SP(idx) = a;
                    b_SP(idx) = b;
                end
                
                % Ring-planet gear contact:
                jdx_RP = find(contains(field_name, 'ring_planet'));
                a_RP = zeros(obj.N_p, 1);
                b_RP = zeros(obj.N_p, 1);
                for idx = 1:obj.N_p
                    jdx = jdx_RP(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    wbull = fitdist(data_idx.load.ov_020.values + eps, 'weibull');
                    
                    a = wbull.a; % scale param.
                    b = wbull.b; % shape param.
                    
                    fprintf('%s:\tRP%d\ta = %e\tb = %e\n', name, idx, a, b);
                    a_RP(idx) = a;
                    b_RP(idx) = b;
                end
                
                a = [a_PL; a_PLC; a_SP; a_RP];
                b = [b_PL; b_PLC; b_SP; b_RP];
                
            else
                % Pinion bearings:
                idx_PIN  = find(contains(field_name, 'pinion_b'));
                n_PIN = length(idx_PIN);
                a_PIN = zeros(size(idx_PIN));
                b_PIN = zeros(size(idx_PIN));
                for idx = 1:n_PIN
                    jdx = idx_PIN(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    
                    [a, b] = obj.bearing(idx).Weibull(data_idx);
                    fprintf('%s\t%s%d\ta = %e\tb = %e\n', name, obj.bearing(idx).name, ...
                                                        idx, a, b);

                    a_PIN(idx) = a;
                    b_PIN(idx) = b;
                end
                
                % Wheel bearings:
                idx_WHE  = find(contains(field_name, 'wheel_b'));
                n_WHE = length(idx_WHE);
                a_WHE = zeros(size(idx_WHE));
                b_WHE = zeros(size(idx_WHE));
                for idx = 1:n_WHE
                    jdx = idx_WHE(idx);
                    data_idx.load = data.load.(field_name{jdx});
                    kdx = n_PIN + idx;
                    
                    [a, b] = obj.bearing(kdx).Weibull(data_idx);
                    fprintf('%s\t%s%d\ta = %e\tb = %e\n', name, obj.bearing(kdx).name, ...
                                                        idx, a, b);

                    a_WHE(idx) = a;
                    b_WHE(idx) = b;
                end
                
                % Pinion-Wheel gear contact:
                jdx_PW = contains(field_name, 'pinion_wheel');
                
                wbull = fitdist(data.load.(field_name{jdx_PW}).ov_020.values + eps, 'weibull');
                
                a = wbull.a; % scale param.
                b = wbull.b; % shape param.
                
                fprintf('%s:\tPW\ta = %e\tb = %e\n', name, a, b);
                
                a = [a_PIN; a_WHE; a];
                b = [b_PIN; b_WHE; b];
            end
            
        end
        
        function save_params(obj, name)
            save(obj.export2struct(), name);
        end
        
    end
    
    methods(Static)
        function [num_cycle, load_edges, speed_mean, N] = LDD(load_signal, load_speed, time_step)
            %LDD produces the load duration distribution of a load_signal
            % with rotational speed load_speed in [1/min.]. time_step is 
            % the inverse of the sampling frequency of both signals.
            % Inputs:
            % - load_signal: load signal to be processed.
            % - load_speed: rotational speed of the element being loaded.
            % - time_step: inverse of the sampling frequency of load and
            % speed signals, given in seconds.
            %
            % Outputs:
            % - num_cycle: cumulative number of cycles
            % - load_edges: edges of the load bins
            % - speed_mean: mean speed at each load bin
            % - N: Number of cycles at each bin
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
            [N, load_edges] = histcounts(load_signal, load_edges);
            
            speed_mean = zeros(1, num_bins - 1);
            for idx = 1:(num_bins - 1)
                if(N(idx) ~= 0.0)
                    range = (load_signal >= load_edges(idx)) & ...
                            (load_signal <= load_edges(idx + 1));
                    speed_mean(idx) = mean(load_speed(range));
                end
            end
            
            % descending order:
            load_edges = fliplr(load_edges);
            speed_mean  = fliplr(speed_mean);
            N          = fliplr(N);
            
            % Assuming that each load is constant during a time step:
            load_time = N*time_step;
            
            % according to [1], (10):
            num_cycle = load_time.*speed_mean/60.0;
            num_cycle = [0.0, cumsum(num_cycle)];
            
            %% Example of plotting:
%             figure;
%             subplot(211)
%             histogram('binEdges', fliplr(load_edges), ...
%                       'binCounts', fliplr(N), ...
%                       'displayStyle', 'stairs', ...
%                       'lineWidth'   , 2.0);
%             title('Histogram');
%             xlabel('Load, [NA]');
%             ylabel('Num. cycles, [-]');
%             
%             subplot(212)
%             histogram('binEdges' , num_cycle, ...
%                       'binCounts', load_edges(2:end), ...
%                       'displayStyle', 'stairs', ...
%                       'lineWidth'   , 2.0);
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
        function obj = set.P_rated(obj, val)
            obj.P_rated = val;
        end
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
        
        function val = get.J_eq(obj)
            % based on: 
            % Borders, J. (2009), Planetary Geartrain Analysis.
            % Accesed on 23.11.2020.
            % http://www.bordersengineering.com/tech_ref/planetary/planetary_analysis.pdf
            %
            
            if(strcmp(obj.configuration, 'parallel'))
                val = obj.J_x(2) + obj.J_x(1)*obj.u^2;
            elseif(strcmp(obj.configuration, 'planetary'))
                r_s = (obj.d(1)*1.0e-3)/2.0;
                r_p = (obj.d(2)*1.0e-3)/2.0;
                r_r = (obj.d(3)*1.0e-3)/2.0;
                aw  =  obj.a_w*1.0e-3;
                
                val = obj.J_x(1) + ...
                      obj.carrier.J_x*(r_s^2)/(2.0*aw*(r_s + r_r)) + ...
                     (obj.J_x(2)/(r_p^2) + obj.mass(2)*aw/(r_s + r_r))*(obj.N_p*r_s^2)/2.0;
            else
                error('prog:input', 'Configuration [%s] is NOT defined.', obj.configuration);
            end
        end
        
        function val = get.k_T(obj)
            % based on: 
            % Borders, J. (2009), Planetary Geartrain Analysis.
            % Accesed on 23.11.2020.
            % http://www.bordersengineering.com/tech_ref/planetary/planetary_analysis.pdf
            %
            
            if(strcmp(obj.configuration, 'parallel'))
                val = -obj.z(2)/obj.z(1);
            elseif(strcmp(obj.configuration, 'planetary'))
                val = obj.d(1)/(4.0*obj.a_w);
            end
        end
    end
end
