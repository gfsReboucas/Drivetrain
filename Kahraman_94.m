classdef Kahraman_94 < Dynamic_Formulation
    %KAHRAMAN_94 Calculates the inertia and stiffness matrices of a
    % multi-stage Drivetrain object according to:
    % [1] A. Kahraman, "Natural Modes of Planetary Gear Trains",
    % Journal of Sound and Vibration, vol. 173, no. 1, pp. 125-130,
    % 1994. https://doi.org/10.1006/jsvi.1994.1222.
    %
    
    methods
        function obj = Kahraman_94(DT)
            obj@Dynamic_Formulation(DT);
        end
        
        %% Calculation:
        function MM = inertia_matrix(obj)
            N = obj.n_DOF(end);
            MM = zeros(N, N);
            
            DT = obj.drive_train;

            MM(1, 1) = DT.J_Rotor;     MM(end, end) = DT.J_Gen;
            
            M_tmp = DT.main_shaft.inertia_matrix("torsional");
            MM(1:2, 1:2) = MM(1:2, 1:2) + M_tmp;
            
            for idx = 1:DT.N_stage
                M_tmp = Kahraman_94.stage_inertia_matrix(DT.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(M_tmp) - 1);
                
                MM(kdx, kdx) = MM(kdx, kdx) + M_tmp;
            end
        end
        
        function KK = stiffness_matrix(obj)
            N = obj.n_DOF(end);
            KK = zeros(N, N);
            
            DT = obj.drive_train;

            K_tmp = DT.main_shaft.stiffness_matrix("torsional");
            KK(1:2, 1:2) = KK(1:2, 1:2) + K_tmp;
            
            for idx = 1:DT.N_stage
                K_tmp = Kahraman_94.stage_stiffness_matrix(DT.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(K_tmp) - 1);
                
                KK(kdx, kdx) = KK(kdx, kdx) + K_tmp;
            end
        end
        
        function DD = damping_matrix(obj)
            N = obj.n_DOF(end);
            DD = zeros(N, N);
            
            DT = obj.drive_train;
            
            D_tmp = DT.main_shaft.damping_matrix("torsional");
            DD(1:2, 1:2) = DD(1:2, 1:2) + D_tmp;
            
            for idx = 1:DT.N_stage
                D_tmp = Kahraman_94.stage_damping_matrix(DT.stage(idx));
                jdx = obj.n_DOF(idx);
                kdx = jdx:(jdx + length(D_tmp) - 1);
                
                DD(kdx, kdx) = DD(kdx, kdx) + D_tmp;
            end
        end
        
        function bb = load_vector(obj)
            bb         = zeros(obj.n_DOF(end), 2);
            bb(1  , 1) = 1.0;
            bb(end, 2) = 1.0;
        end
        
        function nn = calculate_DOF(obj)
            DT = obj.drive_train;
            nn = ones(DT.N_stage + 1, 1)*2.0;
            
            for idx = 1:DT.N_stage
                if(strcmp(DT.stage(idx).configuration, "parallel"))
                    tmp = DT.stage(idx).N_p + 1;
                elseif(strcmp(DT.stage(idx).configuration, "planetary"))
                    tmp = DT.stage(idx).N_p + 2;
                end
                
                nn(idx + 1) = nn(idx) + tmp;
            end
        end
        
        function sol = time_response(obj, varargin)
            %TIME_RESPONSE calculates the time response of the drivetrain.
            % The optional parameters are:
            % - solver: one of MATLAB's ODE solvers. These are recommended
            % for stiff problems: ode15s, ode23s, ode23t, and ode23tb.
            % - time_range: range of integration.
            % - IC: initial conditions.
            % - N_red: Number of modes to be used. The system is solved on
            % (N_red reduced) modal coordinates.
            % - ref_speed: output reference speed, in [1/min].
            % - K_p: Proportional gain.
            % - K_I: Integral gain.
            %
            
            default = {'solver'    , @ode15s, ... 
                       'time_range', [    0.0, ...  % initial time instant, [s]
                                         80.0], ... % final time instant, [s]
                       'IC'        , zeros(2*obj.n_DOF(end), 1), ... % [rad, 1/min]
                       'N_red'     ,      2, ...
                       'ref_speed' ,   1165.9, ... % [1/min]
                       'K_p'       ,   2200.0, ...
                       'K_I'       ,    2200.0};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            solver     = default.solver;
            time_range = default.time_range;
            IC         = default.IC;
            N_red      = default.N_red;
            ref_speed  = default.ref_speed;
            K_p        = default.K_p;
            K_I        = default.K_I;
            
            warning('off', 'Dynamic_Formulation:RB');
            
            % x_red = MS_r * x:
            % Modal reduction, changing coordinates from physical to modal.
            [M_r, K_r, D_r, b_r, Phi_r] = obj.modal_reduction(1:N_red, 3);
            
            % rad/s to 1/min:
            Phi_rv = Phi_r*30.0/pi;
            
            A    = [zeros(N_red), eye(N_red);
                    -M_r\K_r    , -M_r\D_r];
            B_GA = [zeros(N_red, 2);
                    b_r];
            
            B_A = B_GA(:, 1);
            B_G = B_GA(:, 2);
            
            % Back to physical coordinates:
            C = [zeros(1, N_red), Phi_rv(end, :)];
%             C = blkdiag(MS_r*1.0, ... discard angular position later (?)
%                         MS_r*30.0/pi); % and get angular speed in [1/min]
            IC = pinv(blkdiag(Phi_r, Phi_rv))*IC;
            
            % checking controllability using staircase form to avoid
            % problems with ill-conditioned systems.
            [~, ~, ~, ~, k] = ctrbf(A, B_G, C);
            if(sum(k) ~= length(A))
                warning('System is not-controllable.');
            end
            
            % adding integral control:
            A_bar   = [A - K_p*B_G*C, B_G;
                         - K_I    *C, 0.0];
            B_A_bar = [B_A;
                       0.0];
            B_r_bar = [K_p*B_G;
                       K_I];
            IC_bar = [IC;
                      ref_speed];
            
            A_bar   = sparse(A_bar);
            B_A_bar = sparse(B_A_bar);
            B_r_bar = sparse(B_r_bar);
            
            % Smallest natural period:
            [f_n, ~, f_nd, ~, ~] = obj.modal_analysis();
            
            T_n = 1.0/max(f_n);   % undamped
            T_nd = 1.0/max(f_nd); % damped
            time_step = min([T_n, T_nd])/2.0; % initial time-step
            
            opt_ODE = odeset('vectorized' , 'on', ...
                             'jacobian'   , A_bar, ...
                             'initialStep', time_step);

            warning('on', 'Dynamic_Formulation:RB');
            
            data = load('input_load.mat');
            t_aero = data.time;
            T_aero = -data.Mx*1.0e3;
%             T_aero = mean(T_aero);

            T_A = @(t)interp1(t_aero, T_aero, t, 'linear');
%             T_A = @(t) T_aero;
            
            RHS = @(t, x, ref)(A_bar*x + B_A_bar*T_A(t) + B_r_bar*ref_speed);
            
            sol_tmp = solver(@(t, x)RHS(t, x, ref_speed), time_range, IC_bar, opt_ODE);
            
            sol           = sol_tmp;
            sol.t         = sol_tmp.x;
            sol.x_red     = sol_tmp.y(1:end - 1, :);
            sol.x_full    = blkdiag(Phi_r, Phi_rv)*sol.x_red;
            sol.error     = ref_speed - C*sol.x_red;
            sol.T_gen     = K_p*sol.error + sol_tmp.y(end, :);
            sol.ref_speed = ref_speed*ones(size(sol.t));
            sol.T_aero    = T_A(sol.t);
            
            sol = rmfield(sol, {'x', 'y'});
            
%             function x_dot = RHS(t, x)
%                 x_dot = A_bar*x + B_A_bar*T_A(t) + B_r_bar*ref_speed;
%                 x_dot = A*x + B_A*T_A(t) + B_G*TG;
%             end
        end
        
    end
    
    methods(Static)
        function MM = stage_inertia_matrix(stage_idx)
            
            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                MM = zeros(n, n);
                
                m_p = stage_idx.mass(1);
                m_w = stage_idx.mass(2);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                range = 1:2;
                MM(range, range) = diag([m_w*r_w^2 m_p*r_p^2]);
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = stage_idx.N_p + 3;
                MM = zeros(n, n);
                
                m_s = stage_idx.mass(1);
                m_p = stage_idx.mass(2);
                m_c = stage_idx.carrier.mass;
                
                r_c =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                range = 1:(n - 1);
                MM(range, range) = diag([m_c*r_c^2, ...
                                         m_p*r_p^2*ones(1, stage_idx.N_p), ...
                                         m_s*r_s^2]);
            end
            
            range = (n - 1):n;
            
            MM(range, range) = MM(range, range) + ...
                stage_idx.output_shaft.inertia_matrix('torsional');
        end
        
        function KK = stage_stiffness_matrix(stage_idx)

            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                KK = zeros(n, n);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                k_pw = stage_idx.k_mesh;
                
                range = 1:2;
                KK(range, range) = [r_w ^ 2 * k_pw         r_p     * k_pw * r_w;
                                r_w     * k_pw * r_p   r_p ^ 2 * k_pw];
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = stage_idx.N_p + 3;
                KK = zeros(n, n);
                
                % Mesh stiffness:
                % Mesh component:
                sun_pla = stage_idx.sub_set('sun_planet');
                pla_rng = stage_idx.sub_set('planet_ring');

                
                % OBS.: got a negative k_mesh for DTU_10MW planet-ring
                k_sp = sun_pla.k_mesh;
                k_rp = pla_rng.k_mesh;

                r_c =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                KK(    1,     1) =  stage_idx.N_p*(k_rp + k_sp)*r_c^2;
                KK(    1, n - 1) = -stage_idx.N_p*  r_s * k_sp *r_c;
                KK(n - 1,     1) = KK(1, n - 1);
                KK(n - 1, n - 1) =  stage_idx.N_p*k_sp*r_s^2;
                
                for jdx = 2:(n - 2)
                    KK(    1,   jdx) = r_c*r_p*(k_rp - k_sp);
                    KK(jdx  ,     1) = KK(1, jdx);
                    KK(jdx  ,   jdx) = (k_rp + k_sp)*r_p^2;
                    KK(n - 1,   jdx) = r_s*r_p*k_sp;
                    KK(jdx  , n - 1) = KK(n - 1, jdx);
                end
            end
            range = (n - 1):n;
            
            KK(range, range) = KK(range, range) + ...
                stage_idx.output_shaft.stiffness_matrix('torsional');
            
        end
        
        function DD = stage_damping_matrix(stage_idx)

            if(strcmp(stage_idx.configuration, 'parallel'))
                n = 3;
                DD = zeros(n, n);
                
                r_p = (stage_idx.d(1)*1.0e-3)/2.0;
                r_w = (stage_idx.d(2)*1.0e-3)/2.0;
                
                d_pw = 500.0e6;
                
                range = 1:2;
                DD(range, range) = [r_w ^ 2 * d_pw         r_p     * d_pw * r_w;
                                    r_w     * d_pw * r_p   r_p ^ 2 * d_pw];
            elseif(strcmp(stage_idx.configuration, 'planetary'))
                n = stage_idx.N_p + 3;
                DD = zeros(n, n);
                
                d_sp = 500.0e6;
                d_rp = 500.0e6;
                
                r_c =  stage_idx.a_w *1.0e-3;
                r_s = (stage_idx.d(1)*1.0e-3)/2.0;
                r_p = (stage_idx.d(2)*1.0e-3)/2.0;
                
                DD(    1,     1) =  stage_idx.N_p*(d_rp + d_sp)*r_c^2;
                DD(    1, n - 1) = -stage_idx.N_p* r_s  * d_sp *r_c;
                DD(n - 1,     1) = DD(1, n - 1);
                DD(n - 1, n - 1) =  stage_idx.N_p*d_sp*r_s^2;
                
                for jdx = 2:(n - 2)
                    DD(    1,   jdx) = r_c*r_p*(d_rp - d_sp);
                    DD(jdx  ,     1) = DD(1, jdx);
                    DD(jdx  ,   jdx) = (d_rp + d_sp)*r_p^2;
                    DD(n - 1,   jdx) = r_s*r_p*d_sp;
                    DD(jdx  , n - 1) = DD(n - 1, jdx);
                end
            end
            range = (n - 1):n;
            
            DD(range, range) = DD(range, range) + ...
                stage_idx.output_shaft.damping_matrix('torsional');
            
        end
        
    end
end

