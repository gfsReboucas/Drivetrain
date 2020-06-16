classdef Carrier
    %CARRIER
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
    
    properties(Access = private)
        a_w; % [mm], Center distance
        b_g; % [mm], Face width
    end
    
    properties(Dependent)
        d_a;  % [mm],    Tip diameter
        d_f;  % [mm],    Root diameter
        b;    % [mm],    Witdh
        V;    % [mm^3],  Volume
        mass; % [kg],    Mass
        J_x;  % [kg-m^2] Mass moment of inertia (rot. axis)
        J_y;  % [kg-m^2] Mass moment of inertia
        J_z;  % [kg-m^2] Mass moment of inertia
    end
    
    methods
        function obj = Carrier(aw, bg)
            obj.a_w = aw;
            obj.b_g = bg;
        end
        
        function h = plot(obj, varargin)
            if(nargin == 1)
                C = zeros(2, 1);
                plot_prop = {"lineStyle", "-" , "lineWidth", 2.0, "color", [1.0 0.0 0.0]};
            elseif(nargin > 1)
                C = varargin{1};
                plot_prop = varargin(2:end);
            else
                error("prog:input", "Too many variables.");
            end
            
            [x, y] = obj.reference_circle(C);
            
            h = plot(x, y, plot_prop{:});
            
%             x_0 = -[obj.L/2.0 obj.d/2.0]' + C;
%             rectangle('Position', [x_0' obj.L obj.d], plot_prop{:});
%             axis equal;
%             box on;
        end
        
        function h = rectangle(obj, varargin)
            if(nargin == 1)
                C = zeros(2, 1);
                plot_prop = {[1.0 0.0 0.0], "edgeColor", "k", "lineStyle", "-" , "faceColor", [1.0 0.0 0.0]};
            elseif(nargin > 1)
                C = varargin{1};
                plot_prop = varargin(2:end);
            else
                error("prog:input", "Too many variables.");
            end
            
            X = 0.5*obj.b  *[1 -1 -1  1] + C(1);
            Y = 0.5*obj.d_a*[1  1 -1 -1] + C(2);
            h = fill(X, Y, plot_prop{:});

            axis equal;
            box on;
        end
        
        function [x, y] = reference_circle(obj, varargin)
            t = linspace(0.0, 2.0*pi, 1001);
            
            R = obj.a_w;
            
            if(nargin == 1)
                C = zeros(2, 1);
            elseif(nargin == 2)
                C = varargin{1};
            else
                error("prog:input", "Too many variables.");
            end
    
            x = R*cos(t) + C(1);
            y = R*sin(t) + C(2);

            x = [x, C(1)];
            y = [y, C(2)];
            
        end
        
    end
    
    % Get methods:
    methods
        function val = get.d_a(obj)
            val = 2.6*obj.a_w;
        end
        
        function val = get.d_f(obj)
            val = 1.4*obj.a_w;
        end
        
        function val = get.b(obj)
%             val = 0.75*obj.b_g;
            val = 1.2*obj.b_g;
        end
        
        function val = get.V(obj)
            r_out = obj.d_a/2.0;
            r_in  = obj.d_f/2.0;
            area = pi*(r_out^2 - r_in^2);
            val = area*obj.b*1.0e-9;
        end
        
        function val = get.mass(obj)
            rho = Material.rho;
            val = rho*obj.V;
        end
        
        function val = get.J_x(obj)
            r_out = obj.d_a/2.0;
            r_in  = obj.d_f/2.0;
            val = (obj.mass/2.0)*(r_out^2 + r_in^2)*1.0e-6;
        end
        
        function val = get.J_y(obj)
            r_out = obj.d_a/2.0;
            r_in  = obj.d_f/2.0;
            val = (obj.mass/2.0)*(3.0*(r_out^2 + r_in^2) + obj.b^2)*1.0e-6;
        end
        
        function val = get.J_z(obj)
            val = obj.J_y;
        end
    end
    
end

