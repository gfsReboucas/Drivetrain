classdef Dynamic_Formulation
    properties
        drive_train;
        M;     % Inertia matrix
        K;     % Stiffness matrix
        D;     % Damping matix
        n_DOF; % number of degrees of freedom
        load;  % load vector
    end
    
    methods
        function obj = Dynamic_Formulation(DT)
            if(~isa(DT, 'Drivetrain') || ~strcmp(superclasses(DT), 'Drivetrain'))
                error('DT should be a Drivetrain object.');
            end
            
            obj.drive_train = DT;
            
            obj.n_DOF = 2;
            obj.M = obj.inertia_matrix();
            obj.K = obj.stiffness_matrix();
            obj.D = obj.damping_matrix();
            
            obj.load      = ones(obj.n_DOF, 1);
            obj.load(end) = - obj.load(1);
            
        end
        
        function tab = disp(obj)
            
        end
        
        function export_MAT(obj, filename)
            save(filename, 'obj.M', 'obj.K', 'obj.D', 'obj.load');
        end
        
    end
    
    %% Calculation:
    methods
        function [f_n, mode_shape, zeta] = modal_analysis(obj)
            %MODAL_ANALYSIS calculates the resonances and mode shapes of
            % the Drivetrain via a standard algebraic eigenvalue problem [1].
            %
            % [1] D. Inman, Engineering Vibrations, 4th ed. Boston:
            % Pearson, 2014, pp. 402-407.
            %
            
            n = obj.n_DOF(end);
            A = [ zeros(n),     eye(n);
                 -obj.M\obj.K, -obj.M\obj.D];

            [mode_shape, eig_val] = eig(A);
            
            real_part = (eig_val  + eig_val')/2.0;
            imag_part = (eig_val' - eig_val )*sqrt(-1.0)/2;
            
            real_part = diag(real_part); % matrix to vector
            imag_part = diag(imag_part);
            
            omega_n = sqrt(real_part.^2 + imag_part.^2); % lambda to omega_n
            f_n = omega_n./(2.0*pi); % rad/s to Hz
            zeta = -real_part./omega_n;
            
            % sorting in ascending order:
            [f_n, idx] = sort(f_n);
            zeta = zeta(idx);
            mode_shape = mode_shape(:, idx);
            
            % eliminating repeated values:
            f_n = f_n(1:2:end);
            zeta = zeta(1:2:end);
            mode_shape = mode_shape(1:n, 1:2:end);
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(f_n)
                [~, n] = max(abs(mode_shape(:, idx)));

                mode_shape(:, idx) = mode_shape(:, idx)./mode_shape(n, idx);
            end
            
            flag_RB = any(abs(f_n) < 1.0e-2);
            if(flag_RB)
                warning('Dynamic_Formulation:RB', 'Rigid body behavior.');
                idx = (abs(f_n) < 1.0e-2);
                f_n(idx) = 0.0;
                
                f_n = [f_n(~idx);
                       f_n(idx)];
                   
                mode_shape = [mode_shape(:, ~idx), mode_shape(:, idx)];
            end
            
        end
        
        function H = FRF(obj, freq)

            Omega = 2.0*pi*freq;
            i = sqrt(-1.0);
            
            n_Om = length(Omega);
            H = zeros(n_Om, n);
            
            KK = @(x)(obj.K);
            
            for idx = 1:n_Om
                H(idx, :) = (-obj.M.*Omega(idx).^2 + i.*Omega(idx).*obj.D + KK(Omega(idx)))\obj.load;
            end
            
        end
    end
    
    methods(Access = private)
        function MM = inertia_matrix(obj)
            MM = diag([obj.drive_train.J_Rotor;
                       obj.drive_train.J_Gen*obj.drive_train.u^2]);
        end
        
        function KK = stiffness_matrix(obj)
            k_LSS = obj.drive_train.main_shaft.stiffness('torsional');
            k_HSS = obj.drive_train.stage(end).output_shaft.stiffness('torsional');
            
            u2 = obj.drive_train.u^2;
            
            k = (k_LSS*k_HSS*u2)/(k_LSS + k_HSS*u2);
            KK = k*[ 1.0, -1.0;
                    -1.0,  1.0];
        end
        
        function DD = damping_matrix(obj)
            DD = 0.05*obj.stiffness_matrix();
        end
    end
end
