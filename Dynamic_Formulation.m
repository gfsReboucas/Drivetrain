classdef Dynamic_Formulation < Drivetrain
    properties
        M; % Inertia matrix
        K; % Stiffness matrix
        n_DOF; % numbe of degrees of freedom
    end
    
    methods
        function obj = Dynamic_Formulation(DT)
            if(~isa(DT, 'Drivetrain') || ~strcmp(superclasses(DT), 'Drivetrain'))
                error('DT should be a Drivetrain object.');
            end
            
            obj@Drivetrain('N_stage',    DT.N_stage, ...
                           'stage',      DT.stage, ...
                           'n_rotor',    DT.n_rotor, ...
                           'main_shaft', DT.main_shaft, ...
                           'm_Rotor',    DT.m_Rotor, ...
                           'J_Rotor',    DT.J_Rotor, ...
                           'm_Gen',      DT.m_Gen, ...
                           'J_Gen',      DT.J_Gen, ...
                           'S_Hmin',     DT.S_Hmin, ...
                           'S_Fmin',     DT.S_Fmin);
            
            obj.n_DOF = 2;
            obj.M = obj.inertia_matrix();
            obj.K = obj.stiffness_matrix();
        end
        
        function tab = disp(obj)
        end
        
    end
    %% Calculation:
    methods
        function [f_n, mode_shape] = modal_analysis(obj)
            % Cholesky decomposition:
            L = chol(obj.M, "lower");
            K_tilde = L\obj.K/(L');
            
            % correcting numeric erros and make the problem symmetric:
            if(issymmetric(obj.M) && issymmetric(obj.K))
                K_tilde = (K_tilde + K_tilde')/2.0;
            end
            
            [mode_shape, D] = eig(K_tilde);
            D = diag(D);   % matrix to vector
            
            w_n = sqrt(D); % lambda to omega_n

            f_n = w_n./(2.0*pi); % rad/s to Hz
            
            % sorting in ascending order:
            [f_n, idx] = sort(f_n);
            mode_shape = mode_shape(:, idx);
            
            flag_im = any(imag(f_n) ~= 0.0);
            if(flag_im)
                idx = (imag(f_n) ~= 0.0);
                f_n(idx) = 0.0;
                
                f_n = [f_n(~idx);
                       f_n(idx)];
                   
                mode_shape = [mode_shape(:, ~idx), mode_shape(:, idx)];
            end
            
            flag_RB = any(abs(f_n) < 1.0e-2);
            if(flag_RB)
                idx = (abs(f_n) < 1.0e-2);
                f_n(idx) = 0.0;
                
                f_n = [f_n(~idx);
                       f_n(idx)];
                   
                mode_shape = [mode_shape(:, ~idx), mode_shape(:, idx)];
            end
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(f_n)
                [ms_max, n] = max(abs(mode_shape(:, idx)));
                mode_shape(:, idx) = mode_shape(:, idx)*sign(mode_shape(n, idx))./ms_max;
            end
            
        end
    end
    
    methods(Access = private)
        function MM = inertia_matrix(obj)
            MM = diag([obj.J_Rotor;
                       obj.J_Gen*obj.u^2]);
        end
        
        function KK = stiffness_matrix(obj)
            k_LSS = obj.main_shaft.stiffness('torsional');
            k_HSS = obj.stage(end).out_shaft.stiffness('torsional');
            
            uu = obj.u;
            
            k = (k_LSS*k_HSS*uu^2)/(k_LSS + k_HSS*uu^2);
            KK = k*[ 1.0, -1.0;
                    -1.0,  1.0];
        end
    end
end
