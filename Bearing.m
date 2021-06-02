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
        X; % [-], Radial load factor
        Y; % [-], Axial load factor
        a; % [-], Exponent of the life equation
        C; % [N], Basic dynamic load rating
        e; % [-], Force ratio factor
        L_10; % [M], Basic rating life
        
        name(1, :) string;    % [-],       Bearing designation
        type(1, :) string {mustBeMember(type, ["CARB", "CRB", "SRB", "TRB", "none"])} = "none";
                                              % CARB: CARB toroidal roller bearing
                                                      % CRB: Cylindrical Roller Bearing
                                                             % SRB: Spherical Roller Bearing
                                                                    % TRB: Tapered Roller Bearing
        x;       % [mm],      bearing axial position w. r. t. shaft/pin 
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

        D;       % [mm],      Outer diameter
        d;       % [mm],      Inner diameter
        B;       % [mm],      Thickness
    end
    
    methods
        function obj = Bearing(varargin)
            default = {'name'  , 'no_name', ...
                      'type'   , 'none', ...
                      'x'      , 0.0, ...
                      'D'      , 0.0, ...
                      'd'      , 0.0, ...
                      'B'      , 0.0, ...
                      'a'      , 3.0, ...
                      'C'      , 0.0, ...
                      'e'      , 0.0, ...
                      'X'      , 0.0, ...
                      'Y'      , 0.0, ...
                      'L_10'   , 1.0, ...
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
            obj.D       = default.D;
            obj.d       = default.d;
            obj.B       = default.B;
            
            obj.a       = default.a;
            obj.C       = default.C;
            obj.e       = default.e;
            obj.X       = default.X;
            obj.Y       = default.Y;
            obj.L_10    = default.L_10;
            
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
                           "Inner diameter",                   "d",       "mm",        obj.d;
                           "Outer diameter",                   "D",       "mm",        obj.D;
                           "Width",                            "B",       "mm",        obj.B;
                           "Exponent of the life equation",    "a",        "-",        obj.a;
                           "Basic dynamic load rating"    ,    "C",       "N",         obj.C;
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
        
        function data = export2struct(obj)
            warning('off', 'MATLAB:structOnObject');
            data = struct(obj);
            warning('on', 'MATLAB:structOnObject');
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
                              'D'     , obj(1).D, 'd'    , obj(1).d  , 'B'      , obj(1).B);
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
                              'D'      , obj(1).D , 'd'     , obj(1).d   , 'B'      , obj(1).B);
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
        
        function P = dynamic_equiv_load(obj, F_R, F_A)
            cond = F_R./F_A <= obj.e;
            P = (F_R.*obj.X(1) + F_A.*obj.Y(1)).*( cond) + ...
                (F_R.*obj.X(2) + F_A.*obj.Y(2)).*(~cond);
        end
        
        function D = damage_calculation(obj, F_R, F_A, speed, time_step)
            %DAMAGE_CALCULATION performs fatigue damage calculation for a
            % Bearing object. F_R and F_A are the radial and axial forces
            % at the Bearing, which rotates with velocity speed in [1/min.].
            %
            
            % Dynamic equivalent load:
            P = obj.dynamic_equiv_load(F_R, F_A);
            [N, P_edges] = ISO_6336.LDD(P, speed, time_step);
            
            sum_term = N.*power(P_edges, obj.a);
            D = power(obj.L_10*obj.C, -obj.a)*sum(sum_term);
            
        end
        
        function D = damage_analysis(obj, data)
            F_A = data.load.ov_001.values;
            F_y = data.load.ov_002.values;
            F_z = data.load.ov_003.values;
            F_R = sqrt(F_y.^2 + F_z.^2);
            
            % [rad/s] to [1/min.]:
            speed = data.load.ov_017.values.*30.0/pi;
            
            D = obj.damage_calculation(F_R, F_A, speed, data.time_step);
            
        end
        
        function show_LDD(obj, data)
            F_A = data.load.ov_001.values;
            F_y = data.load.ov_002.values;
            F_z = data.load.ov_003.values;
            F_R = sqrt(F_y.^2 + F_z.^2);
            
            % [rad/s] to [1/min.]:
            speed = data.load.ov_017.values.*30.0/pi;
            
            % Dynamic equivalent load:
            P = obj.dynamic_equiv_load(F_R, F_A);
            [N, P_edges] = ISO_6336.LDD(P, speed, data.time_step);
            time_duration = length(speed)*data.time_step;

            histogram('binEdges'    , N*obj.L_10*3.6e3/time_duration, ...
                      'binCounts'   , P_edges(2:end)*1.0e-3, ...
                      'displayStyle', 'stairs', ...
                      'lineWidth'   , 2.0);
        end
        
        function show_histogram(obj, data)
            F_A = data.load.ov_001.values;
            F_y = data.load.ov_002.values;
            F_z = data.load.ov_003.values;
            F_R = sqrt(F_y.^2 + F_z.^2);
            
            % [rad/s] to [1/min.]:
            speed = data.load.ov_017.values.*30.0/pi;
            
            % Dynamic equivalent load:
            P = obj.dynamic_equiv_load(F_R, F_A);
            [N, P_edges, speed_edges] = ISO_6336.LDD(P, speed, data.time_step);
            N = diff(N)*60.0./(speed_edges.*data.time_step);
            N(isnan(N)) = 0;
            
            histogram('binEdges'    , fliplr(P_edges), ...
                      'binCounts'   , fliplr(N), ...
                      'displayStyle', 'stairs', ...
                      'lineWidth'   , 2.0);
        end
        
        function [a, b] = Weibull(obj, data)
            F_A = data.load.ov_001.values;
            F_y = data.load.ov_002.values;
            F_z = data.load.ov_003.values;
            F_R = sqrt(F_y.^2 + F_z.^2);
            
            % Dynamic equivalent load:
            P = obj.dynamic_equiv_load(F_R, F_A);
            
            wbull = fitdist(abs(P) + eps, 'weibull');
            
            a = wbull.a; % scale param.
            b = wbull.b; % shape param.
            
        end
        
    end
    
    methods(Static)
        function data_analysis(tag)
            
            tag = upper(tag);
            switch(tag)
                case 'SRB'
                    type = 'spherical roller bearings';
                    DD_ref =  1030.0;            dd_ref =   710.0;
                    BB_ref =   315.0;            CC_ref = 11164.0;
                    idx = 1;
                case 'CRB1'
                    type = 'cylindrical roller bearing (single row)';
                    DD_ref = 1220.0;            dd_ref = 1000.0;
                    BB_ref =  128.0;            CC_ref = 3690.0;
                    idx = 2;
                case 'CRB2'
                    type = 'cylindrical roller bearing (double row, full complement)';
                    DD_ref =  600.0;            dd_ref =  400.0;
                    BB_ref =  272.0;            CC_ref = 5500.0;
                    idx = 3;
                case 'CARB'
                    type = 'barrel-shaped and toroidal roller bearings';
                    DD_ref = 1030.0;            dd_ref =  710.0;
                    BB_ref =  236.0;            CC_ref = 8800.0;
                    idx = 4;
                case 'TRB1'
                    type = 'tapered roller bearing (single row)';
                    DD_ref =  550.0;            dd_ref =  410.0;
                    BB_ref =   86.0;            CC_ref = 1467.0;
                    idx = 5;
                case 'ALL'
                    tmp = {'SRB', 'CRB1', 'CRB2', 'CARB', 'TRB1'};
                    for idx = 1:numel(tmp)
                        Bearing.data_analysis(tmp{idx});
                    end
                    return;
            end
            
            fig_dim = {'units'   , 'centimeters', ...
                       'position', [5.0 5.0 23.0 16.0]};
            font_style = {'fontName'   , 'Times', ...
                          'fontSize'   , 12.0   , ...
                          'interpreter', 'LaTeX'};
            
            %% Loading data:
            load(sprintf('data\\data_bearing%02d', idx));
            name = {'b', 'C', 'd', 'D', 'Dpw', 'Dw', 'e', 'Lw', 'X1', 'X2', 'Y1', 'Y2', 'Z'};
            BB(BB == 0.0) = nan;            CC(CC == 0.0) = nan;
            dd(dd == 0.0) = nan;            DD(DD == 0.0) = nan;
            Dpw(Dpw == 0.0) = nan;          Dw(Dw == 0.0) = nan;
            Lw(Lw == 0.0) = nan;            ZZ(ZZ == 0.0) = nan;
            data = [BB; CC; dd; DD; Dpw; Dw; ee; Lw; X1; X2; Y1; Y2; ZZ]';
            
%             idx_B   = find(BB);            idx_C   = find(CC);
%             index = intersect(idx_B, idx_C, 'stable');
%             idx_d   = find(dd);            idx_D   = find(DD);
%             index = intersect(index, intersect(idx_d, idx_D, 'stable'), 'stable');
%             idx_Dpw = find(Dpw);           idx_Dw  = find(Dw);
%             index = intersect(index, intersect(idx_Dpw, idx_Dw, 'stable'), 'stable');
%             idx_e   = find(ee);            idx_Lw  = find(Lw);
%             index = intersect(index, intersect(idx_e, idx_Lw, 'stable'), 'stable');
%             idx_X1  = find(X1);            idx_X2  = find(X2);
%             index = intersect(index, intersect(idx_X1, idx_X2, 'stable'), 'stable');
%             idx_Y1  = find(Y1);            idx_Y2  = find(Y2);
%             index = intersect(index, intersect(idx_Y1, idx_Y2, 'stable'), 'stable');
%             idx_ZZ  = find(ZZ);
%             index = intersect(index, idx_ZZ, 'stable');
            
            % Removing empty columns:
            data_var = var(data, 'omitnan');
            data_var(isnan(data_var)) = -1;
            data_var(data_var < eps) = -1;
            idx_0 = find(data_var < 0);
            clear(name{idx_0});
            data(:, idx_0) = [];
            name(:, idx_0) = [];
            
            % Removing Zero elements:
%             index = sort(intersect(intersect(idx_Dw, idx_Dpw, 'stable'), idx_ZZ, 'stable'));
%             data = data(index, :);
            
            % Normalizing:
            data = data./max(data);
            
            % Creating table:
            tab = array2table(data);
            tab.Properties.VariableNames = name;
            
            if(strcmpi(tag, 'CRB1'))
%                 idx_Lw  = find(Lw(index));
%                 data = data(idx_Lw, :);
% %                 
%                 idx_0 = (var(data) < eps);
%                 clear(name{idx_0});
%                 data(:, idx_0) = [];
%                 name(:, idx_0) = [];
%                 
%                 % Updating table:
%                 tab = array2table(data);
%                 tab.Properties.VariableNames = name;
                
                mod_2 = @(b, x)(b(1) + b(2).*x(:, 1));
                fitLw = fitnlm(tab.b, tab.Lw, mod_2, ones(2, 1), 'varNames', {'D', 'L_w'});
                
                figure(fig_dim{:});
                subplot(211)
                hold on;
                scatter(tab.b, tab.Lw, 'filled');
                plot(tab.b, fitLw.feval(tab.b), ...
                     'lineStyle', '-', ...
                     'lineWidth', 2.0);
                xlabel('$B$, [-]'     , font_style{:});
                ylabel('$L_{we}$, [-]', font_style{:});
                
                legend({'data', 'fit'}, 'location', 'best');
                
                subplot(212)
                hold on;
                scatter(tab.b, 35948.0*power(tab.Lw*max(Lw), 8.0/9.0), 'filled');
                scatter(tab.b, 35948.0*power(fitLw.feval(tab.b)*max(Lw), 8.0/9.0), 'filled');
                xlabel('$B$, [-]'           , font_style{:});
                ylabel('$K$, [N/mm$^{8/9}$]', font_style{:});
                title(type);
                
                legend({'K-form.', 'K-fit'}, 'location', 'best');
                
                fig_axes = findobj(gcf, 'Type', 'Axes');
                set(fig_axes, 'fontName', 'Times', 'fontSize', 12.0);
                set(fig_axes, 'box'     , 'on');
            end
            
            %% Correlation analysis:
            figure(fig_dim{:});
            subplot(141)
            plot(tab.d, tab.C, 'ko');
            xlabel('$d$, [-]', font_style{:})
            ylabel('$C$, [-]', font_style{:})
            
            subplot(142)
            plot(tab.D, tab.C, 'ko');
            xlabel('$D$, [-]', font_style{:})
            title(type);
            
            subplot(143)
            plot(tab.b, tab.C, 'ko');
            xlabel('$B$, [-]', font_style{:})
            
            subplot(144)
            plot(tab.Dw, tab.C, 'ko');
            xlabel('$D_w$, [-]', font_style{:})
            
            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'fontName', 'Times', ...
                          'fontSize', 12.0);
            set(fig_axes, 'box'     , 'on');
            
            warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
            figure(fig_dim{:});
            corrplot(tab);
            title(type);
            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'fontName', 'Times', ...
                          'fontSize', 12.0);
            set(fig_axes, 'box'     , 'on');
            
            warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');

            %% Regression analysis:
            x_1ref = DD_ref./max(DD);
            x_2ref = BB_ref./max(BB);
            
            x_1 = tab.D;
            x_2 = tab.b;
            y   = tab.C;
            
            mod_3 = @(b, x)(b(1).*x(:, 2).^2    + ...
                            b(2).*x(:, 1).*x(:, 2));
            b0 = ones(2, 1);
            fitC3 = fitnlm([x_1 x_2], y, mod_3, b0, ...
                           'varNames', {'D', 'B', 'C'});
            fprintf('%s:\n', type);
            fprintf('C_3 ~ %+.3e*B^2 %+.3e*D*B\t with R^2 = %.3e\n', fitC3.Coefficients.Estimate, fitC3.Rsquared.Adjusted);
            fprintf('C_3: (REF.) = %.3e\t(FIT.) = %.3e [N]\n', CC_ref, max(CC)*fitC3.feval(x_1ref, x_2ref));
            
            figure(fig_dim{:});
            hold on;
            scatter3(x_1, x_2, y, 'filled');
            scatter3(x_1, x_2, fitC3.feval(x_1, x_2), 'filled');
            xlabel('$D$, [-]', font_style{:})
            ylabel('$B$, [-]', font_style{:})
            zlabel('$C$, [-]', font_style{:})
            grid on;
            box on;
            set(gca, 'zScale', 'log');
            title(type);
            view(45.0, 45.0);
            
            legend({'data', 'fit'}, 'location', 'best');
            
            fig_axes = findobj(gcf, 'Type', 'Axes');
            set(fig_axes, 'fontName', 'Times', ...
                          'fontSize', 12.0);
            set(fig_axes, 'box'     , 'on');
        end
    end
    
end


