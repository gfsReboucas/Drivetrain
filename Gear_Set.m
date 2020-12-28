classdef Gear_Set < Gear
    % This class implements SOME procedures for the calculation of the load
    % capacity of cylindrical involute gears with external or internal
    % teeth. Specifically the calculation of contact stresses for the
    % assessment of the surface durability of cylindrical gears.
    % In a planetary gear there are two different gear pairs:
    % (1) sun-planet;
    % (2) planet-ring;
    %
    % References:
    % [1] ISO 6336-1:2006 Calculation of load capacity of spur and helical
    % gears -- Part 1: Basic principles, introduction and general influence
    % factors
    % [2] ISO 6336-2:2006 Calculation of load capacity of spur and helical
    % gears -- Part 2: Calculation of surface durability (pitting)
    % [3] ISO/TR 6336-30:2017 Calculation of load capacity of spur and
    % helical gears -- Calculation examples for the application of ISO 6336
    % parts 1, 2, 3, 5
    % [4] Nejad, A. R., Guo, Y., Gao, Z., Moan, T. (2016). Development of a
    % 5 MW reference gearbox for offshore wind turbines. Wind Energy. 
    % https://doi.org/10.1002/we.1884
    % [5] IEC 61400-4:2012 Wind Turbines -- Part 4: Design Requirements for
    % wind turbine gearboxes
    % [6] ISO 21771:2007 Gears -- Cylindrical involute gears and gear pairs
    % -- Concepts and geometry
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
    
    properties(SetAccess = private)
        configuration (1, :) string   {mustBeMember(configuration, ["parallel", "planetary"])} = 'parallel'; % [-], Configuration of the gear set (e.g. parallel, planetary)
        N_p           (1, 1)          {mustBeInteger, mustBePositive}                          = 1;          % [-], Number of planets
        bearing       (1, :) Bearing;                                                                        % [-], Bearing array
    end
    
    properties
        a_w           (1, :)          {mustBeFinite,  mustBePositive}                          = 13;         % [mm], Center distance
        output_shaft  (1, 1) Shaft;                                                                          % [-], Output shaft
    end
    
    properties(Dependent)
        u;             % [-],         Gear ratio
        alpha_wt;      % [deg.],      Working transverse pressure angle
        eps_alpha;     % [-],         Transverse contact ratio
        eps_beta;      % [-],         Overlap ratio
        eps_gamma;     % [-],         Total contact ratio:
        cprime;        % [N/(mm-um)], Maximum single stiffness of a tooth pair
        cprime_th;     % [N/(mm-um)], Theoretical single stiffness
        c_gamma;       % [N/(mm-um)], Mean value of mesh stiffness per unit face witdh
        c_gamma_alpha; % [N/(mm-um)], Mean value of mesh stiffness per unit face witdh (used for K_v, K_Halpha, K_Falpha)
        c_gamma_beta;  % [N/(mm-um)], Mean value of mesh stiffness per unit face witdh (used for K_Hbeta, K_Fbeta)
        k_mesh;        % [N/m],       Mean value of mesh stiffness
        carrier;       % [-],         Planet carrier
        d_Nf;          % [mm],        Start of active profile diameter
        d_Na;          % [mm],        Active tip diameter
        f_pb;          % [um],        Transverse base pitch deviation
        y_alpha;       % [um],        Running-in allowance
        g_alpha;       % [mm],        Length of path of contact
    end
    
    methods
        function obj = Gear_Set(varargin)
            default = {'configuration', 'parallel', ...
                       'z'            , [17   103  ], ...
                       'x'            , [ 0.145 0.0], ...
                       'k'            , [ 1     1  ]*0, ...
                       'bore_ratio'   , [ 1     1  ]*0.5, ...
                       'N_p'          ,   1, ...
                       'a_w'          , 500.0, ...
                       'bearing'      , [Bearing(), Bearing()], ...
                       'shaft'        , Shaft(), ...
                       'm_n'          ,   8.0, ...
                       'alpha_n'      ,  20.0, ...
                       'b'            , 100.0, ...
                       'beta'         ,  15.8, ...
                       'rack_type'    ,  'D', ...
                       'Q'            ,   5.0, ...
                       'R_a'          ,   1.0, ...
                       'material'     , [Material(), Material()]};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            if(length(default.z) < 2)
                error('prog:input', 'There should be at least two gears.');
            elseif(length(default.z) == 3)
                default.z(3) = -abs(default.z(3)); % because the ring is an internal gear
            end
            
            if((length(default.z) ~= length(default.x)) && ...
               (length(default.x) ~= length(default.k)) && ...
               (length(default.k) ~= length(default.bore_ratio)))
                error('prog:input', 'The lengths of z, x, k and bore ratio should be the equal.');
            end
            
            if(std(default.m_n) ~= 0.0)
                error('prog:input', 'Normal modules m_n should be equal for all gears.');
            elseif(std(default.alpha_n) ~= 0.0)
                error('prog:input', 'Pressure angles alpha_n should be equal for all gears.');
            elseif(std(default.beta) ~= 0.0)
                error('prog:input', 'Helix angles beta should be equal for all gears.');
            end
            
            obj@Gear('m_n'       , default.m_n, ...
                     'alpha_n'   , default.alpha_n, ...
                     'type'      , default.rack_type, ...
                     'z'         , default.z, ...
                     'b'         , default.b, ...
                     'x'         , default.x, ...
                     'beta'      , default.beta, ...
                     'k'         , default.k, ...
                     'bore_ratio', default.bore_ratio, ...
                     'Q'         , default.Q, ...
                     'R_a'       , default.R_a, ...
                     'material'  , default.material);
            
            if(strcmp(default.configuration, 'planetary'))
                obj.configuration = default.configuration;
%                 [sun, planet, ring]  =  [1, 2, 3]
                obj.N_p    = default.N_p;
            elseif(strcmp(default.configuration, 'parallel'))
                obj.configuration = default.configuration;
                obj.N_p    = 1;
            else
                error('Gear_Set:configuration_undefined', 'Configuration [%s] is NOT defined.', default.configuration)
            end
            
            obj.a_w = default.a_w;
            obj.bearing = default.bearing;
            obj.output_shaft = default.shaft;
        end
        
        function tab = disp(obj)
            %DISP display some properties of a Gear_Set object

            tmp_vec = NaN(size(obj.z));
            tmp_vec(2) = 1;
            tab_set = {'Gear Ratio',                            'u',       '-',      obj.u  *tmp_vec;
                       'Number of elements'                     'p',       '-',      obj.N_p*tmp_vec;
                       'Normal module',                         'm_n',     'mm',     obj.m_n*tmp_vec;
                       'Normal pressure angle',                 'alpha_n', 'deg.',   obj.alpha_n*tmp_vec;
                       'Helix angle',                           'beta',    'deg.',   obj.beta*tmp_vec;
                       'Face width',                            'b',       'mm',     obj.b*tmp_vec;
                       'Center distance',                       'a_w',     'mm',     obj.a_w*tmp_vec;
                       'Number of teeth',                       'z',       '-',      obj.z;
                       'Profile shift coefficient',             'x',       '-',      obj.x;
                       'Reference diameter',                    'd',       'mm',     obj.d;
                       'Tip diameter',                          'd_a',     'mm',     obj.d_a;
                       'Root diameter',                         'd_f',     'mm',     obj.d_f;
                       'Mass',                                  'm',       'kg',     obj.mass;
                       'Mass moment of inertia (x axis, rot.)', 'J_x',     'kg-m^2', obj.J_x;
                       'Mass moment of inertia (y axis)',       'J_y',     'kg-m^2', obj.J_y;
                       'Mass moment of inertia (z axis)',       'J_z',     'kg-m^2', obj.J_z;
                       'Bearing Names'                          '-+-',     '-',      tmp_vec;
                       'Bearing Types'                          '-+-',     '-',      tmp_vec;
                       };

            Parameter = tab_set(:, 1);
            Symbol    = tab_set(:, 2);
            Unit      = tab_set(:, 3);
            Value     = tab_set(:, 4);
            
            Value = cell2mat(Value);
            Value = num2cell(Value);
            
            Value(cellfun(@isnan,Value)) = {'-+-'}; 
            
            if(strcmp(obj.configuration, 'parallel'))
                for idx = 1:length(obj.z)
                    jdx = 3*idx - 2;
                    Value{end - 1, idx} = join([obj.bearing(jdx:jdx + 2).name], ' / ');
                    Value{end    , idx} = join([obj.bearing(jdx:jdx + 2).type], ' / ');
                end

                v_pinion = Value(:, 1);
                v_wheel  = Value(:, 2);
                
                tab = table(Parameter, Symbol, v_pinion, v_wheel, Unit, ...
                            'variableNames', ["Parameter", "Symbol", "Pinion", "Wheel", "Unit"]);
            elseif(strcmp(obj.configuration, 'planetary'))
                Value{end - 1, 2} = join([obj.bearing(1:2).name], ' / ');
                Value{end    , 2} = join([obj.bearing(1:2).type], ' / ');
                
                Sun     = Value(:, 1);
                Planet  = Value(:, 2);
                Ring    = Value(:, 3);
                Carrier = {'-+-'; '-+-'; '-+-'; ... 
                         '-+-'; '-+-'; '-+-'; ...
                         '-+-'; '-+-'; '-+-'; ...
                         '-+-'; obj.carrier.d_a;
                                obj.carrier.d_f;
                                obj.carrier.mass;
                                obj.carrier.J_x;
                                obj.carrier.J_y;
                                obj.carrier.J_z;
                                join([obj.bearing(3:4).name], ' / ');
                                join([obj.bearing(3:4).type], ' / ')};

                tab = table(Parameter, Symbol, Sun, Planet, Ring, Carrier, Unit);

            else
                error('prog:input', 'Configuration [%s] is NOT defined.', obj.configuration);
            end
            
            if(nargout == 0)
                fprintf('Gear set:\n');
                disp(tab);
                fprintf('Bearings:\n');
                obj.bearing.disp;
                fprintf('Output shaft:\n');
                obj.output_shaft.disp;
                clear tab;
            end
        end
        
        function [tab, tab_str] = comparison(ref, sca)
            tab_str = {'Normal module',                                'm_n',     'mm',     ref.m_n,         sca.m_n,         sca.m_n         / ref.m_n;         % 7
                       'Face width',                                   'b',       'mm',     ref.b,           sca.b,           sca.b           / ref.b;           % 8
                       'Center distance',                              'a_w',     'mm',     ref.a_w,         sca.a_w,         sca.a_w         / ref.a_w;         % 9
                       'Reference diameter (Sun/Pinion)',              'd_1',     'mm',     ref.d(1),        sca.d(1),        sca.d(1)        / ref.d(1);        % 10
                       'Mass (Sun/Pinion)',                            'm_1',     'kg',     ref.mass(1),     sca.mass(1),     sca.mass(1)     / ref.mass(1);     % 11
                       'Mass moment of inertia (Sun/Pinion)',          'J_xx1',   'kg-m^2', ref.J_x(1),      sca.J_x(1),      sca.J_x(1)      / ref.J_x(1);      % 12
                       'Diameter / Output shaft',                      'd',       'mm',     ref.output_shaft.d, sca.output_shaft.d, sca.output_shaft.d / ref.output_shaft.d; % 13
                       'Length / Output shaft',                        'L',       'mm',     ref.output_shaft.L, sca.output_shaft.L, sca.output_shaft.L / ref.output_shaft.L; % 14
                       };
                   
            Parameter = tab_str(:, 1);
            Symbol    = tab_str(:, 2);
            Unit      = tab_str(:, 3);
            Reference = tab_str(:, 4);
            Scale     = tab_str(:, 5);
            Ratio     = tab_str(:, 6);
            
            tab = table(Parameter, Symbol, Scale, Reference, Ratio, Unit);
            
            if(nargout == 0)
                disp(tab);
                clear tab tab_str;
            end
        end
        
        function plot(obj)

            % LINSPECER: Plot lots of lines with very distinguishable and 
            % aesthetically pleasing colors can be dowloaded from MATLAB's
            % File Exchange on:
            % https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap
            color = linspecer(4, 'qualitative');
            
            hold on;
            axis equal;
            box on;
            
            C_p = [obj.a_w, 0.0]';

            if(strcmp(obj.configuration, 'parallel'))
                plot(obj.gear(1), C_p*0.0, 'lineStyle', '-' , 'lineWidth', 2.0, 'color', color(1, :));
                plot(obj.gear(2), C_p*1.0, 'lineStyle', '-' , 'lineWidth', 2.0, 'color', color(2, :));
                
                legend(["Pinion", "Wheel"], 'location', 'best', 'fontName', 'Times', 'fontSize', 12.0);
                
            elseif(strcmp(obj.configuration, 'planetary'))
                plot(obj.gear(1), C_p*0.0, 'lineStyle', '-' , 'lineWidth', 2.0, 'color', color(1, :));
                plot(obj.gear(2), C_p*1.0, 'lineStyle', '-' , 'lineWidth', 2.0, 'color', color(2, :));
                plot(obj.gear(3), C_p*0.0, 'lineStyle', '-' , 'lineWidth', 2.0, 'color', color(3, :));
                plot(obj.carrier, C_p*0.0, 'lineStyle', '-' , 'lineWidth', 2.0, 'color', color(4, :));
                
                RotXY = @(x)[cos(x), sin(x); sin(x), cos(x)];
                
                ang = 2.0*pi/obj.N_p;
                for idx = 2:obj.N_p
                    plot(obj.gear(2), RotXY(ang*(idx - 1))*C_p, 'lineStyle', '-' , 'lineWidth', 2.0, 'color', color(2, :));
                end
                
                legend(["Sun", "Planet", "Ring", "Carrier"], 'location', 'best', 'fontName', 'Times', 'fontSize', 12.0);
            end
            
            xlabel('y');
            ylabel('z');
            
        end
        
        function plot3(obj)
            
            % LINSPECER: Plot lots of lines with very distinguishable and 
            % aesthetically pleasing colors. It can be dowloaded from
            % MATLAB's File Exchange on:
            % https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap
            color = linspecer(4, 'qualitative');
            
            C_p = [obj.a_w, 0.0]';

            if(strcmp(obj.configuration, 'parallel'))
                hp = plot3(obj.gear(1), C_p*0.0, 'edgeColor', 'none', 'lineStyle', '-' , 'faceColor', color(1, :));
                hw = plot3(obj.gear(2), C_p*1.0, 'edgeColor', 'none', 'lineStyle', '-' , 'faceColor', color(2, :));
                
                legend([hp, hw], ["Pinion", "Wheel"], 'location', 'best', 'fontName', 'Times', 'fontSize', 12.0);

            elseif(strcmp(obj.configuration, 'planetary'))
                hs = plot3(obj.gear(1), C_p*0.0, 'edgeColor', 'none', 'lineStyle', '-' , 'faceColor', color(1, :));
                hp = plot3(obj.gear(2), C_p*1.0, 'edgeColor', 'none', 'lineStyle', '-' , 'faceColor', color(2, :));
                hr = plot3(obj.gear(3), C_p*0.0, 'edgeColor', 'none', 'lineStyle', '-' , 'faceColor', color(3, :));
                
                RotXY = @(x)[cos(x), sin(x); sin(x), cos(x)];
                
                ang = 2.0*pi/obj.N_p;
                for idx = 2:obj.N_p
                    plot3(obj.gear(2), RotXY(ang*(idx - 1))*C_p, 'edgeColor', 'none', 'lineStyle', 'none' , 'faceColor', color(2, :));
                end
                
                legend([hs, hp, hr], ["Sun", "Planet", "Ring", "Carrier"], 'location', 'best', 'fontName', 'Times', 'fontSize', 12.0);
            end
            xlabel('y');
            ylabel('x');
            zlabel('z');
            
        end
        
        function data = export2struct(obj)
            warning('off', 'MATLAB:structOnObject');
            data = struct(obj);
            warning('on', 'MATLAB:structOnObject');
            
            data = rmfield(data, 'output_shaft');
            data.output_shaft = obj.output_shaft.export2struct();
            
            if(strcmp(obj.configuration, 'planetary'))
                data = rmfield(data, 'carrier');
                data.carrier = obj.carrier.export2struct();
            end
            
            data = rmfield(data, 'bearing');
            for idx = 1:length(obj.bearing)
                data.bearing(idx) = obj.bearing(idx).export2struct();
            end
            
            data = rmfield(data, 'material');
            for idx = 1:length(obj.material)
                data.material(idx) = obj.material(idx).export2struct();
            end
        end
        
        function rectangle(obj, varargin)
            if(nargin == 1)
                C_0 = zeros(2, 1);
            else
                C_0 = varargin{1};
            end
            
            % LINSPECER: Plot lots of lines with very distinguishable and 
            % aesthetically pleasing colors. It can be dowloaded from
            % MATLAB's File Exchange on:
            % https://se.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap
            color = linspecer(6, 'qualitative');
            
            hold on;
            if(strcmp(obj.configuration, 'parallel'))
                C_w =       [obj.b/2.0 0.0]' + C_0;
                C_p =       [obj.b/2.0 obj.a_w]' + C_0;
                C_s = C_p + [obj.b + obj.output_shaft.L 0.0]'./2.0;
                
                rectangle(obj.gear(1)  , C_p, color(1, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(1, :));
                rectangle(obj.gear(2)  , C_w, color(2, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(2, :));
                rectangle(obj.output_shaft, C_s, color(5, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(5, :));
                
%                 legend([h_p h_w h_s], ['Pinion', 'Wheel', 'Shaft'], 'location', 'best', 'fontName', 'Times', 'fontSize', 12.0);
                
            elseif(strcmp(obj.configuration, 'planetary'))
                C_c = [obj.carrier.b/2.0, 0.0]' + C_0;
                C_p = C_c + [0.0 obj.a_w]';
                C_s = C_c + [obj.carrier.b + obj.output_shaft.L 0.0]'./2.0;

                rectangle(obj.carrier  , C_c, color(4, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(4, :));
                rectangle(obj.gear(3)  , C_c, color(3, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(3, :)); % ring
                rectangle(obj.gear(1)  , C_c, color(1, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(1, :)); % sun
                rectangle(obj.gear(2)  , C_p, color(2, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(2, :)); % planet
                rectangle(obj.output_shaft, C_s, color(5, :), 'edgeColor', 'k', 'lineStyle', '-' , 'faceColor', color(5, :));
                
%                 legend([h_g h_p h_r h_c h_s], ['Sun', 'Planet', 'Ring', 'Carrier', 'Shaft'], 'location', 'best', 'fontName', 'Times', 'fontSize', 12.0);
                
            end
            hold off;
        end
        
        function ks = KISSsoft(obj)
            if(strcmp(obj.configuration, 'parallel'))
                module = 'Z012';
                std_file = 'CylGearPair 1 (spur gear).Z12';
                geo_meth = false; % maybe replace??
            elseif(strcmp(obj.configuration, 'planetary'))
                module = 'Z014';
                std_file = 'PlanetarySet 1 (ISO6336).Z14';
                geo_meth = true;% maybe replace??
            end
            
            ks = KISSsoftCOM(module);
            version = ks.get_version();
            version = strrep(version, '/', '-');            
            
            file_name = sprintf('C:\\Program Files (x86)\\KISSsoft %s\\example\\%s', version, std_file);
            
            try
                ks.load_file(file_name);
            catch err
                delete(ks);
                error(err.identifier, '%s', err.message);
            end
            
            ks.set_var('ZS.AnzahlZwi',         obj.N_p);      % number of planets
            ks.set_var('ZS.Geo.mn'   ,         obj.m_n);      % normal module
            ks.set_var('ZP[0].a'     ,         obj.a_w     ); % center distance
            ks.set_var('ZS.Geo.alfn' , deg2rad(obj.alpha_n)); % normal pressure angle
            ks.set_var('ZS.Geo.beta' , deg2rad(obj.beta)   ); % helix angle
            
            ks.set_var('RechSt.GeometrieMeth', geo_meth);    % tooth geometry according to ISO 21771:2007 % maybe replace??
            
            for idx = 1:numel(obj.z)
                ks.set_var(sprintf('ZR[%d].z'                   , idx - 1), obj.z(idx));
                ks.set_var(sprintf('ZR[%d].x.nul'               , idx - 1), obj.x(idx));
                ks.set_var(sprintf('ZR[%d].b'                   , idx - 1), obj.b     );
                ks.set_var(sprintf('ZR[%d].Tool.type'           , idx - 1), 2);
                ks.set_var(sprintf('ZR[%d].Tool.RefProfile.name', idx - 1), sprintf('1.25 / 0.38 / 1.0 ISO 53:1998 Profil %s', obj.type));
                
                ks.set_var(sprintf('ZR[%d].Vqual', idx - 1), obj.Q);
                ks.set_var(sprintf('ZR[%d].RAH'  , idx - 1), obj.R_a);
                ks.set_var(sprintf('ZR[%d].RAF'  , idx - 1), obj.R_z);
                
                % material properties:
                ks.set_var(sprintf('ZR[%d].mat.bez'    , idx - 1), '18CrNiMo7-6');
                ks.set_var(sprintf('ZR[%d].WerkstArt'  , idx - 1), 'Case-carburized steel');
                ks.set_var(sprintf('ZR[%d].WerkstBh'   , idx - 1), 'case-hardened');
                ks.set_var(sprintf('ZR[%d].mat.comment', idx - 1), 'ISO 6336-5 Figure 9/10 (MQ), Core hardness >=25HRC Jominy J=12mm<HRC28');
            end
            
            if(~ks.calculate())
                delete(ks);
            end
            
            if(nargout == 0)
                delete(ks);
                clear('ks');
            end
        end
        
        function gset = sub_set(obj, option)
            if(~strcmp(obj.configuration, 'planetary'))
                error('Gear_Set:not_planet', 'Not a planetary Gear_Set.');
            end
            
            if(strcmp(option, 'sun_planet'))
                idx = 1:2;
            elseif(strcmp(option, 'planet_ring'))
                idx = 2:3;
            end
            
            gset = Gear_Set('configuration', 'parallel'         , ...
                            'm_n'          , obj.m_n            , ...
                            'alpha_n'      , obj.alpha_n        , ...
                            'z'            , obj.z(idx)         , ...
                            'b'            , obj.b              , ...
                            'beta'         , obj.beta           , ...
                            'x'            , obj.x(idx)         , ...
                            'k'            , obj.k(idx)         , ...
                            'bore_ratio'   , obj.bore_ratio(idx), ...
                            'N_p'          , 1                  , ...
                            'a_w'          , obj.a_w            , ...
                            'rack_type'    , obj.type           , ...
                            'bearing'      , obj.bearing        , ...
                            'shaft'        , obj.output_shaft);
        end
    end
    
    %% Calculation:
    methods
        %% Scaling:
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
            
            m_n_sca = Rack.module(obj_ref.m_n*gamma('m_n'), 'calc', 'nearest');
            
            gamma('m_n') = m_n_sca/obj_ref.m_n;
            
            ref_shaft = obj_ref.output_shaft;
            shaft_sca = Shaft('d'       , ref_shaft.d*gamma('d'), ...
                              'L'       , ref_shaft.L*gamma('L'), ...
                              'bearing' , ref_shaft.bearing, ...
                              'material', ref_shaft.material);
                          
            obj_sca = Gear_Set('configuration', obj_ref.configuration, ...
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
        end
        
        function [obj_sca, gamma, res] = scaled_version(obj_ref, P, n_1, SH_ref, aspect, varargin)
            
            if(strcmp(aspect, 'Gear_Set'))
                n = 4;
            elseif(strcmp(aspect, 'stage'))
                n = 1;
            elseif(strcmp(aspect, 'gear'))
                n = 2;
            else
                error('Aspect = %s is NOT defined.', upper(aspect));
            end
            
            if(~isempty(varargin))
                gamma_0 = varargin{1};
            else
                gamma_0 = ones(n, 1)*0.5;
            end
            
            fun_asp = @(x)(SH_ref - obj_ref.pitting_safety(P, n_1, x, aspect));
            fun_min = @(x)(norm(fun_asp(x))^2);
            
            gamma_min = ones(n, 1)*1.0e-6;
            gamma_Max = ones(n, 1);
            
            fun_ineq = @(x)(1.25 - obj_ref.pitting_safety(P, n_1, x, aspect));

            constraint_fun = @(x)deal(fun_ineq(x), fun_asp(x)); % inequalities, equalities
            
            opt_solver = optimoptions('fmincon', 'display', 'notify');
            
%             id_1 = 'prog:input';
%             id_2 = 'MATLAB:nearlySingularMatrix';
%             warning('off', id_1);
%             warning('off', id_2);
            [gamma, res] = fmincon(fun_min, gamma_0, [], [], [], [], gamma_min, gamma_Max, constraint_fun, opt_solver);
%             warning('on', id_2);
%             warning('on', id_1);
            
            obj_sca = obj_ref.scale_aspect(gamma, aspect);
        end
        
        function [n_new, idx_zero, idx_inp] = set_gear_speed(obj, n_old)
            %SET_GEAR_SPEED sets the speeds of a Gear_Set. The array n_old
            % should specify only the fixed element (if a planetary
            % Gear_Set) and the input element. The other speeds are
            % calculated in this method and should be defined as NaN in
            % n_old.
            %
            
            if(strcmpi(obj.configuration, 'parallel'))
                idx_zero = nan;
                
                U = diag([1.0/obj.u, obj.u]);
                U = flip(U);
            elseif(strcmpi(obj.configuration, 'planetary'))
                idx_zero = find(n_old == 0);
                n_old(idx_zero) = nan;
                zc = num2cell(obj.z);
                [s, p, r] = deal(zc{:});
                r = abs(r);
                U = zeros(4);

                switch(idx_zero)
                    case 1 % fixed: sun
                           % input: ring*
                           % output: carrier* (* or vice-versa)
                        U(:, 4) = [0.0, (r + s)/p, (1.0 + s/r), 1.0]';
                    case 2
                        error('Gear_Set:speed', 'Direct drive not implemented yet.');
                    case 3 % fixed: ring
                           % input: sun*
                           % output: carrier*
                        U(:, 4) = [(1.0 + r/s), (r + s)/p, 0.0, 1.0]';
                    case 4 % fixed: carrier
                           % input: sun*
                           % output: ring*
                        U(:, 3) = [-r/s, r/p, 1.0, 0.0];
                    otherwise
                        if(isempty(idx_zero))
                            error('Gear_Set:speed', 'There are no fixed elements or the size.');
                        elseif(idx_zero > 4)
                            error('Gear_Set:speed', 'There should be a maximum of 4 velocities.');
                        end
                end
                
            end
            
            A = eye(length(U)) - U;
            x = null(A);
            
            idx_inp = find(~isnan(n_old));
            x = x./x(idx_inp);
            n_new = x*n_old(idx_inp);
        end
        
        %% Misc.:
        function k_hat = mesh_stiffness(obj, varargin)
            %MESH_STIFFNESS Calculates the mesh stiffness for a Gear_Set
            % according to [1-2].
            % [1] Gu, X., Velex, P., Sainsot, P., and Bruyère, J. (June 1,
            % 2015). "Analytical Investigations on the Mesh Stiffness
            % Function of Solid Spur and Helical Gears." ASME. J. Mech.
            % Des. June 2015; 137(6): 063301. 
            % https://doi.org/10.1115/1.4030272
            % [2] Velex Philippe (April 11th 2012). On the Modelling of 
            % Spur and Helical Gear Dynamic Behaviour, Mechanical 
            % Engineering, Murat Gokcek, IntechOpen, DOI: 10.5772/36157.
            %
            
            default = {'t',       0.0, ...
                       'k',       5, ...
                       'Omega_1', 1.0};
            
            default = scaling_factor.process_varargin(default, varargin);
            t       = default.t;
            k       = default.k;
            Omega_1 = default.Omega_1;
            
            % Meshing period, [2]:
            T_m = pi*obj.m_n*cosd(obj.alpha_n)/(2.0*obj.d_b(1)*Omega_1);
            tau = t/T_m; % normalized time instant
            
            eps_a = obj.eps_alpha;            eps_b = obj.eps_beta;
            
            % Eq. (10), [1]:
            Xi = @(k)((0.7 + 0.09./(k.*eps_a).^2).*sinc(k.*eps_a) + ...
                - cos(pi.*k.*eps_a).*0.09./(k.*eps_a).^2);
            
            gen_term = @(kk, tt)(Xi(kk).*sinc(kk.*eps_b).* ...
                                 cos(pi.*kk.*(2.0.*tt - eps_a - eps_b)));
            
            k_hat = 1.0;
            for idx = 1:k
                k_hat = k_hat + 2.0*gen_term(idx, tau);
            end
            
            k_hat = k_hat.*obj.k_mesh;
            
        end
        
        function m_tot = get_mass(obj)
            if(strcmp(obj.configuration, 'parallel'))
                m_tot = sum(obj.mass);
            elseif(strcmp(obj.configuration, 'planetary'))
                m_tot = sum([1.0, obj.N_p, 1.0].*obj.mass) + obj.carrier.mass;
            end
        end
        
        function J_z = get_J_z(obj)
            % with respect to gear 1.
            if(strcmp(obj.configuration, 'parallel'))
                J_z = sum(obj.J_z) + obj.mass(2)*obj.a_w^2;
            elseif(strcmp(obj.configuration, 'planetary'))
                J_z = sum(obj.J_z) + obj.N_p*obj.mass(2)*obj.a_w^2 + obj.carrier.J_z;
            end
        end
        
        function aw = find_center_distance(obj, alpha_wt_star)
            aw = abs(obj.z(1) + obj.z(2))*obj.m_n*cosd(obj.alpha_t)/(2.0*cosd(alpha_wt_star)*cosd(obj.beta));
        end
        
        function kinematic_tree(obj)
            % +--------------+----------------------------------+--------------------------------+-----------------------+
            % |              |              Stage 1             |             Stage 2            |        Stage 3        |
            % |    Element   +------------+---------------------+------------+-------------------+------------+----------+
            % |              |    Prev.   |         Next        |    Prev.   |        Next       |    Prev.   |   Next   |
            % +--------------+------------+---------------------+------------+-------------------+------------+----------+
            % | Output shaft |    frame   | sun   1 / carrier 2 |    frame   | sun   2 / Wheel 3 |    frame   | Pinion 3 |
            % +--------------+------------+---------------------+------------+-------------------+------------+----------+
            % |  Sun/Pinion  | O. shaft 1 |          NA         | O. shaft 2 |         NA        | O. shaft 3 |    NA    |
            % +--------------+------------+---------------------+------------+-------------------+------------+----------+
            % | Planet/Wheel |    pin 1   |          NA         |    pin 2   |         NA        | O. shaft 2 |    NA    |
            % +--------------+------------+---------------------+------------+-------------------+------------+----------+
            % |     Ring     |    Frame   |          NA         |    Frame   |         NA        |                       |
            % +--------------+------------+---------------------+------------+-------------------+                       |
            % |    Carrier   |  M. shaft  |        pin 1        | O. shaft 1 |       pin 2       |                       |
            % +--------------+------------+---------------------+------------+-------------------+                       |
            % |      pin     |  carrier 1 |       planet 1      |  carrier 2 |      planet 2     |                       |
            % +--------------+------------+---------------------+------------+-------------------+-----------------------+
            
            clr = ['r', 'g', 'b', 'k', 'c']';
            
            origin = 's';
            input  = 'o';
            output = '>';
            
            if(strcmp(obj.configuration, 'parallel'))
                
            elseif(strcmp(obj.configuration, 'planetary'))
                SUN = struct;
                PLT = struct;
                RNG = struct;
                CAR = struct; % 2
                SFT = struct; % 1
                HSG = struct;
                
                HSG.y = 0.0;
                
                
                SUN.y = 5.0;
                SUN.out = obj.carrier.b/2.0;
                
                PLT.y = 4.0;
                PLT.A = -0.77*obj.b;
                PLT.B = -PLT.A;
                
                RNG.y = 3.0;
                
                SFT.y   = 1.0;
                SFT.inp = -obj.output_shaft.L/2.0;
                SFT.out = -SFT.inp;
                
                CAR.y = 2.0;
                CAR.A   = -0.6*obj.carrier.b/2.0;
                CAR.B   =  0.6*obj.carrier.b/2.0;
                CAR.inp =  CAR.A;
                CAR.out =  0.0;
                
                figure;
                subplot(121)
                hold on;
                plot(0.0    , SUN.y, [clr(1, :), origin]);
                plot(SUN.out, SUN.y, [clr(1, :), output]);
                
                plot(0.0    , CAR.y, [clr(4, :), origin])
                plot(CAR.inp, CAR.y, [clr(4, :), input ]);
                
                plot(0.0    , SFT.y, [clr(5, :), origin]);
                plot(SFT.inp, SFT.y, [clr(5, :), input ]);
                plot(SFT.out, SFT.y, [clr(5, :), output]);
                
                set(gca, 'ytick'     ,  0:5);
                set(gca, 'yticklabel', ["HSG", "SFT", "CAR", "RNG", "PLT", "SUN"]);
                set(gca, 'ylim', [0 5]);
                
                subplot(122)
                hold on;
                
            end
        end
    end
    
    %% Set methods:
    methods
        function obj = set.a_w(obj, val)
            obj.a_w = val;
        end
        
        function obj = set.cprime(obj, val)
            obj.cprime = val;
        end
        
        function obj = set.output_shaft(obj, val)
            obj.output_shaft = val;
        end
        
    end
    
    %% Get methods:
    methods
        function val = get.k_mesh(obj)
            val = obj.c_gamma.*obj.b.*(1.0e6);
        end
        
        function val = get.carrier(obj)
            if(strcmp(obj.configuration, 'planetary'))
                val = Carrier(obj.a_w, obj.b);
            else
                error('prog:input', 'Only Planetary gear sets have a planet carrier.');
            end
        end
        
        function g = gear(obj, idx)
            g = Gear('m_n'        , obj.m_n,  ...
                     'alpha_n'    , obj.alpha_n,  ...
                     'type'       , obj.type,  ...
                     'z'          , obj.z(idx),   ...
                     'b'          , obj.b, ...
                     'x'          , obj.x(idx),  ...
                     'beta'       , obj.beta,  ...
                     'k'          , obj.k(idx),  ...
                     'bore_ratio' , obj.bore_ratio(idx),  ...
                     'Q'          , obj.Q,  ...
                     'R_a'        , obj.R_a);
        end
        
        function val = get.u(obj)
            if(strcmp(obj.configuration, 'parallel'))
                val = obj.z(2)/obj.z(1);
            elseif(strcmp(obj.configuration, 'planetary'))
                val = 1.0 + abs(obj.z(3))/obj.z(1);
            else
                error('prog:input', 'Configuration [%s] is NOT defined.', obj.configuration);
            end
        end
        
        function val = get.alpha_wt(obj)
            num = obj.m_n*cosd(obj.alpha_t);
            den = 2.0*obj.a_w*cosd(obj.beta);
            val = acosd(abs(obj.z(1) + obj.z(2))*num/den);
        end
        
        function val = get.eps_alpha(obj)
            xi_Nfw1 = zeros(3, 1);
            xi_Nfw2 = zeros(3, 1);
            
            % roll angles from the root form diameter to the working pitch point, 
            % limited by the:
            % (1) base diameters: Eq. (33)
            xi_Nfw1(1) = tand(obj.alpha_wt);
            xi_Nfw2(1) = xi_Nfw1(1);
            
            % (2) root form diameters: Eq. (34-35)
            xi_Nfw1(2) = xi_Nfw1(1) - tan(acos(obj.d_b(1)/obj.d_Nf(1)));
            xi_Nfw2(2) = xi_Nfw2(1) - tan(acos(obj.d_b(2)/obj.d_Nf(2)));
            
            % (3) active tip diameters of the wheel/pinion: Eq. (36-37)
            xi_Nfw1(3) = (tan(acos(obj.d_b(2)/obj.d_Na(2))) - xi_Nfw1(1))*obj.z(2)/obj.z(1);
            xi_Nfw2(3) = (tan(acos(obj.d_b(1)/obj.d_Na(1))) - xi_Nfw2(1))*obj.z(1)/obj.z(2);
            
            xi_Nfw1(xi_Nfw1 < 0.0) = [];
            xi_Nfw2(xi_Nfw2 < 0.0) = [];

            xi_Nfw1 = min(xi_Nfw1);
            xi_Nfw2 = min(xi_Nfw2);
            
            if(isempty(xi_Nfw1) || isempty(xi_Nfw2))
                error('Auxiliary coefficients xi_Nfw are all lower than 0.');
            end

            % roll angle from the working pitch point to the active tip diameter: Eq. (38)
            z21 = abs(obj.z(2))/obj.z(1); % to account for the sign of z in ring gears.
            xi_Naw1 = xi_Nfw2*z21;

            % pinion angular pitch:
            tau_1 = 2.0*pi/obj.z(1);

            % Eq. (32)
            val = (xi_Nfw1 + xi_Naw1)/tau_1;
        end
        
        function val = get.eps_beta(obj)
            bb = mean(obj.b);
            val = bb*sind(obj.beta)/(pi*obj.m_n);
        end
        
        function val = get.eps_gamma(obj)
            val = obj.eps_alpha + obj.eps_beta;
        end
        
        function val = get.cprime(obj)
            % Correction factor:
            C_M = 0.8; % solid disk gears

            % Gear blank factor:
            C_R = 1.0; % solid disk gears

            % Basic rack factor:
            alpha_Pn = obj.alpha_n; % [rad.], Normal pressure angle of basic rack
            C_B1 = (1.0 + 0.5*(1.2 - obj.h_fP/obj.m_n))*(1.0 - 0.02*(20.0 - alpha_Pn)); % 0.975
            C_B2 = (1.0 + 0.5*(1.2 - obj.h_fP/obj.m_n))*(1.0 - 0.02*(20.0 - alpha_Pn));

            C_B = 0.5*(C_B1 + C_B2); % 0.975
            
            val = obj.cprime_th*C_M*C_R*C_B*cosd(obj.beta);
        end
        
        function val = get.cprime_th(obj)
            C_1 =  0.04723;      C_2 =  0.15551;      C_3 =  0.25791;
            C_4 = -0.00635;      C_5 = -0.11654;      C_6 = -0.00193;
            C_7 = -0.24188;      C_8 =  0.00529;      C_9 =  0.00182;
            
            % q' is the minimum value for the flexibility of a pair of teeth
            qprime = C_1 + C_2/obj.z_n(1) + C_3/obj.z_n(2) + C_4*obj.x(1) + C_5*(obj.x(1)/obj.z_n(1)) + ...
                C_6*obj.x(2) + C_7*(obj.x(2)/obj.z_n(2)) + C_8*obj.x(1)^2 + C_9*obj.x(1)^2; % [mm-um/N]
            
            % c'_th is the theoretical single stiffness:
            val = 1.0/qprime;

        end
        
        function val = get.c_gamma(obj)
            val = obj.c_gamma_alpha + obj.c_gamma_beta;
        end
        
        function val = get.c_gamma_alpha(obj)
            val = obj.cprime*(0.75*obj.eps_alpha + 0.25);
        end
        
        function val = get.c_gamma_beta(obj)
            val = 0.85*obj.c_gamma_alpha;
        end
        
        function val = get.d_Nf(obj)
            % ISO 21771, Sec. 5.4.1, Eqs. (64-67)
            
            d_Fa = obj.d_a;
            val(1) = sqrt((2.0*obj.a_w*sind(obj.alpha_wt) - sign(obj.z(2))*sqrt(d_Fa(2)^2 - obj.d_b(2)^2))^2 + obj.d_b(1)^2);
            val(2) = sqrt((2.0*obj.a_w*sind(obj.alpha_wt) -                sqrt(d_Fa(1)^2 - obj.d_b(1)^2))^2 + obj.d_b(2)^2);
            
            if(strcmp(obj.configuration, 'planetary'))
                val(3) = nan;
            end
            
            for idx = 1:length(obj.z)
                if(obj.d_Ff(idx) > val(idx))
                    val(idx) = obj.d_Ff(idx);
                end
            end
        end
        
        function val = get.d_Na(obj)
            % ISO 21771, Sec. 5.4.1, Eqs. (68-69)
            
            aw = 2.0*obj.a_w*sind(obj.alpha_wt);
            if(obj.d_Nf(1) == obj.d_Ff(1))
                dNa2 = sqrt((aw -                sqrt(obj.d_Ff(1)^2 - obj.d_b(1)^2))^2 + obj.d_b(2)^2);
            else
                dNa2 = obj.d_a(2); % d_Fa
            end
            
            if(obj.d_Nf(2) == obj.d_Ff(2))
                dNa1 = sqrt((aw - sign(obj.z(2))*sqrt(obj.d_Ff(2)^2 - obj.d_b(2)^2))^2 + obj.d_b(1)^2);
            else
                dNa1 = obj.d_a(1); % d_Fa
            end
            
            val = [dNa1, dNa2];
            if(strcmp(obj.configuration, 'planetary'))
                val(3) = nan;
            end

        end
        
        function val = get.f_pb(obj)
            val = max(obj.f_pt.*cosd(obj.alpha_t));
        end
        
        function val = get.y_alpha(obj)
            if(obj.f_pb >= 40.0) % [um]
                val = 3.0; % [um]
            else
                val = obj.f_pb*75.0e-3;
            end
        end
        
        function val = get.g_alpha(obj)
            val = (sqrt(obj.d_Na(1)^2 - obj.d_b(1)^2) + sign(obj.z(2))*(sqrt(obj.d_Na(1)^2 - obj.d_b(1)^2) - 2.0*obj.a_w*sind(obj.alpha_wt)))/2.0;
        end
        
    end
    
    %% Validation:
    methods(Static)
        function test_k_mesh()
            % to do: test obj.mesh_stiffness against data from [1], Fig. 5.
            %
            % [1] Gu, X., Velex, P., Sainsot, P., and Bruyère, J. (June 1,
            % 2015). "Analytical Investigations on the Mesh Stiffness
            % Function of Solid Spur and Helical Gears." ASME. J. Mech.
            % Des. June 2015; 137(6): 063301. 
            % https://doi.org/10.1115/1.4030272
            % [2] Rohatgi, A. (July, 2020). WebPlotDigitizer
            % https://automeris.io/WebPlotDigitizer
            %
            % +-----------------------+--------------------+
            % |          tau          |   k(tau), Eq. (9)  |
            % +-----------------------+--------------------+
            % | 5.551115123125783e-17 | 0.9171063829787234 |
            % +-----------------------+--------------------+
            % |  0.04885993485342027  | 0.9337872340425531 |
            % +-----------------------+--------------------+
            % |  0.10097719869706845  | 0.9516595744680851 |
            % +-----------------------+--------------------+
            % |  0.15146579804560267  | 0.9698723404255318 |
            % +-----------------------+--------------------+
            % |   0.2035830618892509  | 0.9885957446808511 |
            % +-----------------------+--------------------+
            % |   0.252442996742671   | 1.0078297872340425 |
            % +-----------------------+--------------------+
            % |  0.30293159609120535  |  1.027744680851064 |
            % +-----------------------+--------------------+
            % |  0.35342019543973946  | 1.0394893617021277 |
            % +-----------------------+--------------------+
            % |   0.4039087947882737  | 1.0420425531914894 |
            % +-----------------------+--------------------+
            % |   0.4543973941368078  |  1.043404255319149 |
            % +-----------------------+--------------------+
            % |   0.504885993485342   |  1.044936170212766 |
            % +-----------------------+--------------------+
            % |   0.5553745928338762  | 1.0456170212765958 |
            % +-----------------------+--------------------+
            % |   0.6058631921824105  | 1.0454468085106383 |
            % +-----------------------+--------------------+
            % |   0.6563517915309447  | 1.0442553191489363 |
            % +-----------------------+--------------------+
            % |   0.7068403908794789  | 1.0296170212765958 |
            % +-----------------------+--------------------+
            % |   0.7557003257328991  |  1.009531914893617 |
            % +-----------------------+--------------------+
            % |   0.8078175895765471  | 0.9892765957446809 |
            % +-----------------------+--------------------+
            % |   0.8599348534201954  | 0.9693617021276595 |
            % +-----------------------+--------------------+
            % |   0.9104234527687296  | 0.9506382978723404 |
            % +-----------------------+--------------------+
            % |   0.9609120521172638  | 0.9314042553191488 |
            % +-----------------------+--------------------+
            % |   1.011400651465798   | 0.9201702127659575 |
            % +-----------------------+--------------------+
            % |   1.0602605863192183  |  0.937531914893617 |
            % +-----------------------+--------------------+
            % |   1.1123778501628665  | 0.9552340425531914 |
            % +-----------------------+--------------------+
            % |   1.1644951140065147  | 0.9739574468085106 |
            % +-----------------------+--------------------+
            % |   1.211726384364821   | 0.9928510638297873 |
            % +-----------------------+--------------------+
            % |   1.2638436482084692  | 1.0119148936170212 |
            % +-----------------------+--------------------+
            % |   1.3159609120521174  | 1.0316595744680852 |
            % +-----------------------+--------------------+
            % |   1.3631921824104234  | 1.0398297872340425 |
            % +-----------------------+--------------------+
            % |   1.4153094462540716  | 1.0423829787234042 |
            % +-----------------------+--------------------+
            % |   1.4641693811074918  |  1.043744680851064 |
            % +-----------------------+--------------------+
            % |   1.514657980456026   |  1.045276595744681 |
            % +-----------------------+--------------------+
            % |   1.5667752442996743  | 1.0456170212765958 |
            % +-----------------------+--------------------+
            % |   1.6172638436482085  | 1.0451063829787235 |
            % +-----------------------+--------------------+
            % |   1.6677524429967427  | 1.0442553191489363 |
            % +-----------------------+--------------------+
            % |   1.7182410423452772  | 1.0253617021276595 |
            % +-----------------------+--------------------+
            % |   1.7687296416938108  |  1.005276595744681 |
            % +-----------------------+--------------------+
            % |   1.8192182410423452  | 0.9851914893617021 |
            % +-----------------------+--------------------+
            % |   1.8713355048859934  | 0.9656170212765958 |
            % +-----------------------+--------------------+
            % |   1.9218241042345279  | 0.9472340425531914 |
            % +-----------------------+--------------------+
            % |   1.970684039087948   | 0.9274893617021276 |
            % +-----------------------+--------------------+
            %
            
%             data = [5.551115123125783e-17	0.9171063829787234
%                     0.04885993485342027     0.9337872340425531
%                     0.10097719869706845     0.9516595744680851
%                     0.15146579804560267     0.9698723404255318
%                     0.2035830618892509      0.9885957446808511
%                     0.252442996742671       1.0078297872340425
%                     0.30293159609120535     1.027744680851064
%                     0.35342019543973946     1.0394893617021277
%                     0.4039087947882737      1.0420425531914894
%                     0.4543973941368078      1.043404255319149
%                     0.504885993485342       1.044936170212766
%                     0.5553745928338762      1.0456170212765958
%                     0.6058631921824105      1.0454468085106383
%                     0.6563517915309447      1.0442553191489363
%                     0.7068403908794789      1.0296170212765958
%                     0.7557003257328991      1.009531914893617
%                     0.8078175895765471      0.9892765957446809
%                     0.8599348534201954      0.9693617021276595
%                     0.9104234527687296      0.9506382978723404
%                     0.9609120521172638      0.9314042553191488
%                     1.011400651465798       0.9201702127659575
%                     1.0602605863192183      0.937531914893617
%                     1.1123778501628665      0.9552340425531914
%                     1.1644951140065147      0.9739574468085106
%                     1.211726384364821       0.9928510638297873
%                     1.2638436482084692      1.0119148936170212
%                     1.3159609120521174      1.0316595744680852
%                     1.3631921824104234      1.0398297872340425
%                     1.4153094462540716      1.0423829787234042
%                     1.4641693811074918      1.043744680851064
%                     1.514657980456026       1.045276595744681
%                     1.5667752442996743      1.0456170212765958
%                     1.6172638436482085  	1.0451063829787235
%                     1.6677524429967427      1.0442553191489363
%                     1.7182410423452772      1.0253617021276595
%                     1.7687296416938108      1.005276595744681
%                     1.8192182410423452      0.9851914893617021
%                     1.8713355048859934      0.9656170212765958
%                     1.9218241042345279      0.9472340425531914
%                     1.970684039087948       0.9274893617021276];
            
            % Reading data from [1], Fig. 5. The data was extracted using
            % WebPlotDigitizer [2].
            mat = readmatrix('Fig_5_Gu_Velex_2015.csv');
            tau9  = mat(:, 1);          tau9  = tau9(1:41);
            eq_9  = mat(:, 2);          eq_9  = eq_9(1:41);
            t_LDP = mat(:, 3);          LDP   = mat(:, 4);
            t_VSA = mat(:, 5);          VSA   = mat(:, 6);          
            tau10 = mat(:, 7);          tau10 = tau10(1:41);
            eq_10 = mat(:, 8);          eq_10 = eq_10(1:41);
            
            figure();
            subplot(121)
            hold on;
            plot(tau9, eq_9 ,  'marker', 'x' , 'markerSize', 10.0, 'lineStyle', 'none', 'markerFaceColor',  [0.346 0.536 0.691]);
            plot(tau10, eq_10, 'marker', 'o' , 'markerSize',  5.0, 'lineStyle', 'none', 'markerFaceColor',  [0.915 0.281 0.287]);
            plot(t_VSA, VSA, 'lineStyle', ':' , 'lineWidth',  2.0, 'color', [0.441 0.749 0.432]);
            plot(t_LDP, LDP, 'lineStyle', '-.', 'lineWidth',  2.0, 'color', [1.0   0.598 0.2  ]);
            xlim([0.0 2.0]);
            ylim([0.9 1.1]);
            
            xlabel('\tau');
            ylabel('k(\tau)');
            legend({'Eq. (9)', 'Eq. (10)', 'VSA', 'LDP'}, 'location', 'best')
            set(gca, 'box', 'on');
            
        end
    end
    
end