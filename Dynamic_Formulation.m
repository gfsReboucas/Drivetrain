classdef Dynamic_Formulation
    %DYNAMIC_FORMULATION Generates the dynamic matrices for a Drivetrain 
    % object, based on its Gear_Set and inertia elements. The equations of
    % motion are given as:
    % M q'' + (D + Omega_c G) q' + (K_m + K_b - Omega_c^2 K_Omega) q =
    %       = Omega_c^2 c + b [T, F],
    % where:
    % - q = q(t) is the (angular/translational) displacement vector, whose
    % first and second time derivatives are q' and q'', respectively.
    % - M is the inertia matrix.
    % - D is the damping matrix.
    % - G is the gyroscopic matrix.
    % - K_m is the mesh stiffness matrix.
    % - K_b is the bearing stiffness matrix.
    % - K_Omega is the centripetal stiffness matrix.
    % - c is the centripetal force vector.
    % - b is the external load vector.
    %
   
    properties
        drive_train;
        n_DOF;   % number of degrees of freedom
        M;       % Inertia matrix
        D;       % Damping matrix
        G;       % Gyroscopic matrix
        K_m;     % Mesh stiffness matrix
        K_b;     % Bearing stiffness matrix
        K_Omega; % Centripetal stiffness matrix
        c;       % Centripetal force vector
        b;       % External load vector
        A;       % State matrix
        DOF_description; 
    end
    
    methods
        function obj = Dynamic_Formulation(DT)
            if(~isa(DT, 'Drivetrain') || ~strcmp(superclasses(DT), 'Drivetrain'))
                error('DT should be a Drivetrain object.');
            end
            
            obj.drive_train = DT;
            
            obj.n_DOF = obj.calculate_DOF();
            
            [obj.M, ...
             obj.G]       = obj.inertia_matrix();
            obj.D         = obj.damping_matrix();
            [obj.K_m, ...
             obj.K_b, ...
             obj.K_Omega] = obj.stiffness_matrix();
            obj.c         = obj.centripetal_force_vector();
            obj.b         = obj.external_load_vector();
            obj.A         = obj.state_matrix();
            
            obj.DOF_description = obj.explain_DOF();
        end
        
        function tab = disp(obj)
            
        end
        
        function desc = explain_DOF(obj)
            desc = cell(obj.n_DOF(end), 1);
            desc{1} = 'Rotor angular displacement';
            desc{2} = 'Generator angular displacement';
        end
    end
    
    %% Calculation:
    methods
        function [f_n, MS, f_nd, MS_d, zeta] = modal_analysis(obj)
            %MODAL_ANALYSIS calculates the resonances and mode shapes of
            % the Drivetrain via a standard algebraic eigenvalue problem [1].
            %
            % [1] D. Inman, Engineering Vibrations, 4th ed. Boston:
            % Pearson, 2014, pp. 402-407.
            %
            
            %% Undamped analysis:
            L = chol(obj.M, 'lower');
            K_tilde = L\(obj.K_m + obj.K_b)/(L');
            [MS, EV] = eig(K_tilde);
            MS = (L')\MS;
            
            EV = diag(EV); % matrix to vector
            omega_n = sqrt(EV); % lambda to omega
            
            % sorting in ascending order:
            [omega_n, idx] = sort(omega_n);
            MS = MS(:, idx);
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(omega_n)
                [~, n] = max(abs(MS(:, idx)));

                MS(:, idx) = MS(:, idx)./MS(n, idx);
            end
            
            %% Damped case:
            [MS_d, EV_d] = eig(obj.A);
            
            real_part = (EV_d  + EV_d')/2.0;
            imag_part = (EV_d' - EV_d )*sqrt(-1.0)/2;
            
            real_part = diag(real_part); % matrix to vector
            imag_part = diag(imag_part);
            
            omega_n = sqrt(real_part.^2 + imag_part.^2); % lambda to omega_n
            f_n = omega_n./(2.0*pi); % rad/s to Hz
            zeta = -real_part./omega_n;
            
            % sorting in ascending order:
            [f_n, idx] = sort(f_n);
            zeta = zeta(idx);
            MS_d = MS_d(:, idx);
            
            % eliminating repeated values:
            f_n = f_n(1:2:end);
            zeta = zeta(1:2:end);
            f_nd = f_n.*sqrt(1.0 - zeta.^2);
            n = obj.n_DOF(end);
            MS_d = MS_d(1:n, 1:2:end);
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(f_n)
                [~, n] = max(abs(MS_d(:, idx)));

                MS_d(:, idx) = MS_d(:, idx)./MS_d(n, idx);
            end
            
            flag_RB = any(abs(f_n) < 1.0e-2);
            if(flag_RB)
                warning('Dynamic_Formulation:RB', 'Rigid body behavior.');
                idx = (abs(f_n) < 1.0e-2);
                f_n(idx) = 0.0;
                
                f_n = [f_n(~idx);
                        f_n(idx)];
                   
                MS_d = [MS_d(:, ~idx), MS_d(:, idx)];
            end
            
        end
        
        function [M_r, D_r, K_mr, K_br, b_r, MS_r] = modal_reduction(obj, N, prec)
            [~, MS] = obj.modal_analysis();
            [M_r, D_r, K_mr, K_br, b_r, MS_r] = ...
                Dynamic_Formulation.modal_truncation(obj.M, obj.D, obj.K_m, obj.K_b, obj.b, MS, N, prec);
            
        end
        
        function sol = time_response(obj, varargin)
            %TIME_RESPONSE calculates the time response of the drivetrain.
            % The optional parameters are:
            % - solver: one of MATLAB's ODE solvers. These are recommended
            % for stiff problems: ode15s, ode23s, ode23t, and ode23tb.
            % - time_range: range of integration.
            % - IC: initial conditions.
            % - ref_speed: output reference speed, in [1/min].
            % - K_p: Proportional gain.
            % - K_I: Integral gain.
            % - N_modes: Number of modes to be used. The system is solved
            % using modal coordinates.
            %
            
            IC = zeros(2*obj.n_DOF(end), 1);
            
            default = {'solver'    , @ode23t  , ... 
                       'time_range', [    0.0 , ...  % initial time instant, [s]
                                        100.0], ... % final time instant, [s]
                       'IC'        ,    IC    , ... % [rad, 1/min]
                       'K_p'       ,   2200.0 , ...
                       'K_I'       ,    220.0 , ...
                       'f_sample'  ,    200.0};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            solver     = default.solver;
            time_range = default.time_range;
            IC         = default.IC;
            K_p        = default.K_p;
            K_I        = default.K_I;
            T_sample   = 1.0/default.f_sample;
            
            time = time_range(1):T_sample:time_range(end);
            
            warning('off', 'Dynamic_Formulation:RB');
            
            MM = obj.M;
            DD = obj.D;
            Km = obj.K_m;
            Kb = obj.K_b;
            bb = obj.b;
            
            KK = Km + Kb;
            
            N_modes = obj.n_DOF(end);
            I_n     = eye(N_modes);

            % rad/s to 1/min:
            I_nv = I_n*30.0/pi;
            
            AA = blkdiag(MM, MM);
            BB = [zeros(N_modes),  MM;
                        -KK     , -DD];
            b_GA = [zeros(N_modes, 2);
                        bb];
                
            b_A = b_GA(:, 1);
            b_G = b_GA(:, 2);
            
            % To get the generator speed:
            C = [zeros(1, N_modes), I_nv(end, :)];

            % checking controllability using staircase form to avoid
            % problems with ill-conditioned systems.
            [~, ~, ~, ~, k] = ctrbf(AA\BB, b_G, C);
            if(sum(k) ~= length(AA))
                warning('System is not-controllable.');
            end
            
            % coupling PI control:
            A_bar = blkdiag(AA, 1.0);
            B_bar = [BB - K_p*b_G*C, b_G;
                        - K_I    *C, 0.0];
            c_A_bar = [b_A;
                       0.0];
            c_r_bar = [K_p*b_G;
                       K_I];

            A_bar   = sparse(A_bar);
            B_bar   = sparse(B_bar);
            c_A_bar = sparse(c_A_bar);
            c_r_bar = sparse(c_r_bar);
            
            % Smallest natural period:
            [f_n, ~, f_nd] = obj.modal_analysis();
            
            T_n = 1.0/max(f_n);   % undamped
            T_nd = 1.0/max(f_nd); % damped
            time_step = min([T_n, T_nd])/2.0; % initial time-step
            
            opt_ODE = odeset('vectorized' , 'on', ...
                             'mass'       , A_bar, ...
                             'jacobian'   , B_bar, ...
                             'initialStep', time_step, ...
                             'relTol'     , 1.0e-6, ...
                             'absTol'     , 1.0e-8);
%                              'maxStep'    , 1.0e-3, ...

            warning('on', 'Dynamic_Formulation:RB');
            
            data = load('input_load.mat');
            t_load = data.time;
            T_aero = data.Mx*1.0e3;
            gen_speed = data.w_GEN*(30.0/pi); % from [rad/s] to [1/min]
%             k = find(t_load >= 0.2*time_range(end), 1, 'first');
%             transient = [linspace(0.0, 1.0, k), ...
%                          ones(1, length(t_load) - k)]';
%             T_aero = T_aero.*transient;
            T_A = @(t)interp1(t_load, T_aero, t, 'linear');
            ref_gen_speed = @(t)interp1(t_load, gen_speed, t, 'linear');

            IC_bar = [IC;
                      0.0];

            RHS = @(t, x)(B_bar*x + c_A_bar*T_A(t) + c_r_bar*ref_gen_speed(t));
            
            [~, x_solver] = solver(@(t, x)RHS(t, x), time, IC_bar, opt_ODE);
            x_solver = x_solver';
            nn = obj.n_DOF(end);
            range = 1:nn;
            
            vel_acc = A_bar\RHS(time, x_solver);
            
            sol           = struct;
            sol.solver    = solver; 
            sol.t         = time;
            sol.x_solver  = x_solver; 
            sol.x         = blkdiag(I_n, I_nv)*sol.x_solver(1:end - 1, :);
            sol.gen_speed = C*sol.x_solver(1:end - 1, :);
            sol.ref_speed = ref_gen_speed(sol.t);
            sol.error     = sol.ref_speed - sol.gen_speed;
            sol.T_gen     = K_p*sol.error + x_solver(end, :);
            sol.T_aero    = T_A(sol.t);
            sol.M         = MM;
            sol.D         = DD;
            sol.K         = KK;
            sol.b         = bb;
            sol.K_p       = K_p;
            sol.K_I       = K_I;
            sol.pos       = sol.x  (range     , :);
            sol.vel       = sol.x  (range + nn, :);
            sol.acc       = vel_acc(range + nn, :);
%             sol.load      = sol.K*sol.pos + ...
%                             sol.D*sol.vel + ...
%                             sol.M*sol.acc;
            sol.DOF_description = obj.DOF_description;
            
%             function x_dot = RHS(t, x)
%                 %   A_bar *
%                 x_dot = B_bar*x + B_A_bar*T_A(t) + B_r_bar*ref_gen_speed(t);
%             end
        end
        
        function sol = Newmark_integration(obj,varargin)
            
            x_0 = zeros(obj.n_DOF(end), 1);
            
            default = {'time_range', [    0.0 , ...  % initial time instant, [s]
                                        100.0], ... % final time instant, [s]
                       'x_0'       ,    x_0   , ... % initial position
                       'v_0'       ,    x_0   , ... % initial velocity
                       'K_p'       ,   2200.0 , ...
                       'K_I'       ,    220.0 , ...
                       'f_sample'  ,    200.0 , ...
                       'option'    , 'average'};
            
            default = scaling_factor.process_varargin(default, varargin);
            
            time_range = default.time_range;
            x_0        = default.x_0;
            v_0        = default.v_0;
            K_p        = default.K_p;
            K_I        = default.K_I;
            T_sample   = 1.0/default.f_sample;
            option     = default.option;
            
            delta = 1.0/2.0;
            if(strcmpi(option, 'average')) % unconditionally stable
                alpha = 1.0/4.0;
            elseif(strcmpi(option, 'linear'))
                alpha = 1.0/6.0;
            end
            
            time = time_range(1):T_sample:time_range(2);
            
            n = length(x_0);     nt = length(time);
            
            data = load('input_load.mat');
            t_load = data.time;
            T_aero = data.Mx*1.0e3;
            gen_speed = data.w_GEN*(30.0/pi); % from [rad/s] to [1/min]

            T_A = interp1(t_load, T_aero, time, 'linear');
            ref_gen_speed = interp1(t_load, gen_speed, time, 'linear');

            warning('off', 'Dynamic_Formulation:RB');
            
            C = zeros(1, n);
            C(end) = 30.0/pi;
            
            output = zeros(1, nt);
            error = zeros(1, nt);
            T_G = zeros(1, nt);
            
            pos = zeros(n, nt);
            vel = zeros(n, nt);
            acc = zeros(n, nt);
            
            a_1 = 1.0/(alpha*T_sample^2);
            a_2 = delta/(alpha*T_sample);
            a_3 = 1.0/(alpha*T_sample);
            a_4 = 1.0/(2.0*alpha) - 1.0;
            a_5 = delta/alpha - 1.0;
            a_6 = (T_sample/2.0)*(delta/alpha - 2.0);
            a_7 = T_sample*(1 - delta);
            a_8 = delta*T_sample;
            
            KK = obj.K_b + obj.K_m;
            K_hat = a_1*obj.M + a_2*obj.D + KK;
            
            %% Initial calculation:
            idx = 1;

            pos(:, idx) = x_0;
            vel(:, idx) = v_0;
            
            output(:, idx) = C*vel(:, idx);
            error(:, idx) = ref_gen_speed(:, idx) - output(:, idx);
            sum_error = error(:, idx);
            T_G(:, idx) = K_p*error(:, idx) + K_I*T_sample*sum_error;
            T_AG = [T_A(:, idx);
                    T_G(:, idx)];

            acc(:, idx) = obj.M\(obj.b*T_AG - obj.D*vel(:, idx) - KK*pos(:, idx));
            
            %% for each time step:
            nt = nt - 1;
            
            for idx = 1:nt
                Delta_R_hat = obj.b*T_AG + obj.M*(a_1*pos(:, idx) + a_3*vel(:, idx) + a_4*acc(:, idx)) + ...
                                           obj.D*(a_2*pos(:, idx) + a_5*vel(:, idx) + a_6*acc(:, idx));
                
                pos(:, idx + 1) = K_hat\Delta_R_hat;
                acc(:, idx + 1) = a_1*(pos(:, idx + 1) - pos(:, idx)) - a_3*vel(:, idx) - a_4*acc(:, idx);
                vel(:, idx + 1) = vel(:, idx) + a_7*acc(:, idx) + a_8*acc(:, idx + 1);
                
                output(:, idx + 1) = C*vel(:, idx + 1);
                error (:, idx + 1) = ref_gen_speed(:, idx + 1) - output(:, idx + 1);
                sum_error = sum_error + error(:, idx + 1);
                T_G(:, idx + 1) = K_p*error(:, idx + 1) + K_I*T_sample*sum_error;
                T_AG = [T_A(:, idx + 1);
                        T_G(:, idx + 1)];
                            
            end
            
            sol = struct;
            sol.solver = 'Dynamic_Formulation.Newmark';
            sol.delta  = delta;
            sol.alpha  = alpha;
            sol.t      = time;
            sol.pos    = pos;
            sol.vel    = vel*C(end);
            sol.acc    = acc;
            sol.M      = obj.M;
            sol.D      = obj.D;
            sol.K      = KK;
            sol.C      = C;
            sol.b      = obj.b;
            sol.K_p    = K_p;
            sol.K_I    = K_I;

        end
        
        function [Hx, Hv, Ha] = FRF(obj, freq)
            % FRF Calculates the Frequency Response Function of the model
            % for the frequencies freq in [Hz], returns:
            % - Hx: Receptance function
            % - Hv: Mobility function
            % - Ha: Accelerance function
            %

            Omega = 2.0*pi*freq;
            i = sqrt(-1.0);
            
            n_Om = length(Omega);
            Hx = zeros([n_Om, size(obj.b)]);
            Hv = Hx;        Ha = Hx;
            KK = obj.K_m + obj.K_b;
            
            for idx = 1:n_Om
                Hx(idx, :, :) = (-obj.M.*Omega(idx).^2 + i.*Omega(idx).*obj.D + KK)\obj.b;
                Hv(idx, :, :) = Hx(idx, :, :).*Omega(idx);
                Ha(idx, :, :) = Hv(idx, :, :).*Omega(idx);
            end
            
        end
        
        function [MM, GG] = inertia_matrix(obj)
            MM = diag([obj.drive_train.J_Rotor;
                       obj.drive_train.J_Gen*obj.drive_train.u^2]);
            GG = zeros(obj.n_DOF(end));
        end
        
        function DD = damping_matrix(obj)
            [~, Kb, ~] = obj.stiffness_matrix();
            DD = 0.05*Kb;
        end
        
        function [Km, Kb, KOm] = stiffness_matrix(obj)
            Km = zeros(obj.n_DOF(end));
            KOm = zeros(obj.n_DOF(end));
            
            % Bearing matrix:
            k_LSS = obj.drive_train.main_shaft.stiffness('torsional');
            k_HSS = obj.drive_train.stage(end).output_shaft.stiffness('torsional');
            
            u2 = obj.drive_train.u^2;
            
            k = (k_LSS*k_HSS*u2)/(k_LSS + k_HSS*u2);
            Kb = k*[ 1.0, -1.0;
                    -1.0,  1.0];
        end
        
        function cc = centripetal_force_vector(obj)
            cc = zeros(obj.n_DOF(end), 1);
        end
        
        function bb = external_load_vector(obj)
            bb = eye(obj.n_DOF(end));
        end
        
        function nn = calculate_DOF(obj)
            nn = 2;
        end
        
        function AA = state_matrix(obj)
            KK = obj.K_m + obj.K_b;
            
            n = obj.n_DOF(end);
            
            AA = [ zeros(n),  eye(n);
                  -obj.M\KK, -obj.M\obj.D];
        end
    end
    
    methods(Static)
        function [M_r, D_r, K_mr, K_br, b_r, Phi_r] = modal_truncation(MM, DD, Km, Kb, bb, Phi, N, prec)
            %MODAL_TRUNCATION 
            % u = Phi_r * u_r
            %
            
            Phi_r = Phi(:, N);
            
            eps = power(10.0, -prec);
            
            M_r = Phi_r'*MM*Phi_r;    M_r( abs(M_r)  < eps) = 0.0;
            D_r = Phi_r'*DD*Phi_r;    D_r( abs(D_r)  < eps) = 0.0;
            K_mr = Phi_r'*Km*Phi_r;   K_mr(abs(K_mr) < eps) = 0.0;
            K_br = Phi_r'*Kb*Phi_r;   K_br(abs(K_br) < eps) = 0.0;
            b_r = Phi_r'*bb;          b_r( abs(b_r)  < eps) = 0.0;
            
        end
                
        function [sol, time, pos, vel, acc] = Newmark(time, x0, v0, M, D, K, R, option)
            %NEWMARK solution of linear mechanical systems using Newmark's
            % method according to [1] Table 9.3.
            %
            % Parameters:
            % time: time interval for integration
            % x0: initial position
            % v0: initial velocity
            % M: Inertia matrix
            % D: Damping matrix
            % K: Stiffness matrix
            % R: vector of externally applied loads
            % option: average or linear version of the method
            %
            % [1] K. J. Bathe, Finite Element Procedures, 2nd Edition, 2014.
            % Avaliable online at: (access 24.09.2020)
            % http://web.mit.edu/kjb/www/Books/FEP_2nd_Edition_4th_Printing.pdf
            %
            
            delta = 1.0/2.0;
            if(strcmpi(option, 'average')) % unconditionally stable
                alpha = 1.0/4.0;
            elseif(strcmpi(option, 'linear'))
                alpha = 1.0/6.0;
            end
            
            Delta_t = time(2) - time(1);
            n = length(x0);     nt = length(time);
            
            pos = zeros(n, nt);
            vel = zeros(n, nt);
            acc = zeros(n, nt);
            
            a_1 = 1.0/(alpha*Delta_t^2);
            a_2 = delta/(alpha*Delta_t);
            a_3 = 1.0/(alpha*Delta_t);
            a_4 = 1.0/(2.0*alpha) - 1.0;
            a_5 = delta/alpha - 1.0;
            a_6 = (Delta_t/2.0)*(delta/alpha - 2.0);
            a_7 = Delta_t*(1 - delta);
            a_8 = delta*Delta_t;
            
            K_hat = a_1*M + a_2*D + K;
            
            %% Initial calculation:
            idx = 1;
            pos(:, idx) = x0;
            vel(:, idx) = v0;
            acc(:, idx) = M\(R(:, idx) - D*vel(:, idx) - K*pos(:, idx));
            
            %% for each time step:
            nt = nt - 1;
            
            for idx = 1:nt
                Delta_R_hat = R(:, idx + 1) + M*(a_1*pos(:, idx) + a_3*vel(:, idx) + a_4*acc(:, idx)) + ...
                                              D*(a_2*pos(:, idx) + a_5*vel(:, idx) + a_6*acc(:, idx));
                
                pos(:, idx + 1) = K_hat\Delta_R_hat;
                acc(:, idx + 1) = a_1*(pos(:, idx + 1) - pos(:, idx)) - a_3*vel(:, idx) - a_4*acc(:, idx);
                vel(:, idx + 1) = vel(:, idx) + a_7*acc(:, idx) + a_8*acc(:, idx + 1);
            end
            
            sol = struct;
            sol.solver = 'Newmark';
            sol.delta = delta;
            sol.alpha = alpha;
            sol.t = time;
            sol.x = pos;
            sol.v = vel;
            sol.a = acc;

        end
        
        function [sol, time, pos, vel, acc] = Wilson(time, x0, v0, M, D, K, p)
            
            theta = 1.42; % optimal accuracy
                        
            Ts = time(2) - time(1);
            n = length(x0);     nt = length(time);
            
            pos = zeros(n, nt);
            vel = zeros(n, nt);
            acc = zeros(n, nt);
            
            Delta_p = diff(p, [], 2);

            r_1 = 6.0/(theta*Ts)^2;
            r_2 = 3.0/(theta*Ts);
            r_3 = 6.0/(theta*Ts);
            r_4 = theta*Ts/2.0;
            r_5 = Ts/2.0;
            r_6 = (Ts^2)/2.0;
            r_7 = (Ts^2)/6.0;
            
            K_hat = r_1*M + r_2*D + K;
            a     = r_3*M + 3.0*D;
            b     = 3.0*M + r_4*D;
            
            %% Initial calculation:
            idx = 1;
            pos(:, idx) = x0;
            vel(:, idx) = v0;
            acc(:, idx) = M\(p(:, idx) - D*vel(:, idx) - K*pos(:, idx));
            
            %% for each time step:
            nt = nt - 1;
            
            for idx = 1:nt
                delta_p_hat = theta*Delta_p(:, idx) + a*vel(:, idx) + b*acc(:, idx);
                
                delta_pos = K_hat\delta_p_hat;
                delta_acc = r_1*delta_pos - r_3*vel(:, idx) - 3.0*acc(:, idx);
                
                Delta_acc = delta_acc*(1.0/theta);
                Delta_vel = Ts*acc(:, idx) + r_5*Delta_acc;
                Delta_pos = Ts*vel(:, idx) + r_6*acc(:, idx) + r_7*Delta_acc;
                
                pos(:, idx + 1) = pos(:, idx) + Delta_pos;
                vel(:, idx + 1) = vel(:, idx) + Delta_vel;
                acc(:, idx + 1) = acc(:, idx) + Delta_acc;
            end
            
            sol = struct;
            sol.solver = 'Wilson';
            sol.theta = theta;
            sol.t = time;
            sol.x = pos;
            sol.v = vel;
            sol.a = acc;

        end
        
        function [sol, time, pos, vel, acc] = Bathe(time, x0, v0, M, D, K, R)
            %BATHE solution of linear mechanical systems using Bathe's
            % method according to [1] Table 9.4.
            %
            % Parameters:
            % time: time interval for integration
            % x0: initial position
            % v0: initial velocity
            % M: Inertia matrix
            % D: Damping matrix
            % K: Stiffness matrix
            % R: vector of externally applied loads
            %
            % [1] K. J. Bathe, Finite Element Procedures, 2nd Edition, 2014.
            % Avaliable online at: (access 24.09.2020)
            % http://web.mit.edu/kjb/www/Books/FEP_2nd_Edition_4th_Printing.pdf
            %
            
            Delta_t = time(2) - time(1);
            n = length(x0);     nt = length(time);
            
            pos = zeros(n, nt);
            vel = zeros(n, nt);
            acc = zeros(n, nt);
            
            a_1 = 16.0/(Delta_t^2);
            a_2 =  4.0/ Delta_t;
            a_3 =  9.0/(Delta_t^2);
            a_4 =  3.0/ Delta_t;
            a_5 =  2.0*a_2;
            a_6 = 12.0/(Delta_t^2);
            a_7 = -3.0/(Delta_t^2);
            a_8 = -1.0/ Delta_t;
            
            K_hat1 = K + a_1*M + a_2*D;
            K_hat2 = K + a_3*M + a_4*D;
            
            R_mid = (R(:, 1:end-1) + R(:, 2:end))/2.0;
            
            %% Initial calculation:
            idx = 1;
            pos(:, idx) = x0;
            vel(:, idx) = v0;
            acc(:, idx) = M\(R(:, idx) - D*vel(:, idx) - K*pos(:, idx));
            
            %% for each time step:
            nt = nt - 1;
            
            for idx = 1:nt
                % First sub-step:
                R_hat1 = R_mid(:, idx) + M*(a_1*pos(:, idx) + a_5*vel(:, idx) + acc(:, idx)) + ...
                                         D*(a_2*pos(:, idx) +     vel(:, idx));
                pos_1 = K_hat1\R_hat1;
                vel_1 = a_2*(pos_1 - pos(:, idx)) - vel(:, idx);
%                 acc_1 = a_2*(vel_1 - vel(:, idx)) - acc(:, idx);
                
                % Second sub-step:
                R_hat2 = R(:, idx + 1) + M*(a_6*pos_1 + a_7*pos(:, idx) + a_2*vel_1 + a_8*vel(:, idx)) + ...
                                         D*(a_2*pos_1 + a_8*pos(:, idx));
                pos(:, idx + 1) = K_hat2\R_hat2;
                vel(:, idx + 1) = -a_8*pos(:, idx) - a_2*pos_1 + a_4*pos(:, idx + 1);
                acc(:, idx + 1) = -a_8*vel(:, idx) - a_2*vel_1 + a_4*vel(:, idx + 1);
            end
            
            sol = struct;
            sol.solver = 'Bathe';
            sol.t = time;
            sol.x = pos;
            sol.v = vel;
            sol.a = acc;

        end
        
        function test_Newmark()
            
            time = 0.0:0.1:1.0;
            x0 = 0.0;
            v0 = 0.0;
            M = 0.2533;
            K = 10.0;
            D = 0.1592;
            p = zeros(size(time));
            
            p(2:6) = [5.0, 8.6602, 10.0, 8.6603, 5.0];
            option = 'linear';
            [~, time, pos, vel, acc] = Dynamic_Formulation.Newmark(time, x0, v0, M, D, K, p, option);
            
            table(time', acc', vel', pos', 'variableNames', ["t_i", "a_i", "v_i", "x_i"])
        end
        
        function test_Newmark_02()
            ex_93 = [0.0 0.00673 0.0505 0.189 0.485 0.961 1.58 2.23 2.76 3.00 2.85 2.28 1.40;
                     0.0 0.364   1.35   2.68  4.00  4.95  5.34 5.13 4.48 3.64 2.90 2.44 2.31];
            
            n = 2;
            M = diag([2.0, 1.0]);

            K = diag([6.0, 4.0]);
            K(1, 2) = -2.0;     K(2, 1) = K(1, 2);
            
            D = zeros(n);
            
            r = [0.0, 10.0]';
            
            Ts = 0.28;
            t = 0.0:Ts:(12*Ts);
            x0 = zeros(n, 1);
            v0 = x0;
            
            R = repmat(r, 1, length(t));
            
            option = 'average';
            [~, time, pos, ~, ~] = Dynamic_Formulation.Newmark(t, x0, v0, M, D, K, R, option);
            rel_diff = 100*(ex_93'-pos')./ex_93';
            tab = [(0:12)', ...
                   time', ...
                   ex_93(1, :)', ... 
                   pos(1, :)', ... 
                   rel_diff(:, 1), ...
                   ex_93(2, :)', ... 
                   pos(2, :)', ...
                   rel_diff(:, 2)];

            array2table(tab, 'variableNames', ["idx", "t", "U_1(ref)", "U_1(num)", "rel_diff_1_pc", "U_2(ref)", "U_2(num)", "rel_diff_2_pc"])

        end
        
        function test_Wilson()
            n = 5; % degrees of freedom
            m = 208.6; % [kip-s^2/ft]
            EI = 5.469e10; % [kip-ft^2]
            h = 120.0; % [ft]
            p = 1.0e3; % [kip]
            
            time = 0.0:0.1:2.0;
            x0 = zeros(n, 1);
            v0 = x0;

            M = m*diag([1.0 1.0 1.0 1.0 0.5]);
            
            K_d = diag([18.83 14.65 14.06 9.878 1.608]);
            K_u = zeros(n);
            K_u(1, 2) = -11.90;     K_u(1, 3) =   4.773;    K_u(1, 4) = -1.193;     K_u(1, 5) =  0.1989;
                                    K_u(2, 3) = -10.71;     K_u(2, 4) =  4.177;     K_u(2, 5) = -0.6961;
                                                            K_u(3, 4) = -9.514;     K_u(3, 5) =  2.586;
                                                                                    K_u(4, 5) = -3.646;
            K = K_u + K_u' + K_d;
            K = K*EI/(h^3);

            D = zeros(n);
            
            P = zeros(n, 1);
            P(n) = p;
            P = repmat(P, 1, length(time));
            P(n, 1) = 0.0;
            
            [~, time, pos, ~, ~] = Dynamic_Formulation.Wilson(time, x0, v0, M, D, K, P);
            
            array2table([time', pos'], 'variableNames', ["t_i", "x_1", "x_2", "x_3", "x_4", "x_5"])
            
        end
        
        function test_Bathe()
            %TEST_BATHE replicates the numerical results presented in [1],
            % Example 9.4 and 9.7 (analytical).
            %
            % see also: Bathe, test_Newmark_02
            %
            ex_94 = [0.0 0.00458 0.0445 0.183 0.486 0.979 1.62 2.28 2.81 3.03 2.83 2.21 1.28;
                     0.0 0.373   1.38   2.73  4.04  4.97  5.31 5.06 4.38 3.55 2.85 2.46 2.40];
            n = 2;
            M = diag([2.0, 1.0]);

            K = diag([6.0, 4.0]);
            K(1, 2) = -2.0;     K(2, 1) = K(1, 2);
            
            % Analytical solution according to Example 9.7, [1]:
%             MS = [1.0,  1.0/2.0;
%                   1.0, -1.0]*diag([1.0/sqrt(3.0), sqrt(2.0/3.0)]);
%             ana_sol = @(t) (MS*[(5.0/sqrt(3.0)    )*( 1.0 - cos(t*sqrt(2.0)));
%                                 (2.0*sqrt(2.0/3.0))*(-1.0 + cos(t*sqrt(5.0)))]);
            
            D = zeros(n);
            
            r = [0.0, 10.0]';
            
            Ts = 0.28;
            t = 0.0:Ts:(12*Ts);
            x0 = zeros(n, 1);
            v0 = x0;
            
            R = repmat(r, 1, length(t));
            
            [~, time, pos, ~, ~] = Dynamic_Formulation.Bathe(t, x0, v0, M, D, K, R);
            
            rel_diff = 100*(ex_94'-pos')./ex_94';
            tab = [(0:12)', ...
                   time', ...
                   ex_94(1, :)', ... 
                   pos(1, :)', ... 
                   rel_diff(:, 1), ...
                   ex_94(2, :)', ... 
                   pos(2, :)', ...
                   rel_diff(:, 2)];

            array2table(tab, 'variableNames', ["idx", "t", "U_1(ref)", "U_1(num)", "rel_diff_1_pc", "U_2(ref)", "U_2(num)", "rel_diff_2_pc"])
            
        end
        
    end
    
end
