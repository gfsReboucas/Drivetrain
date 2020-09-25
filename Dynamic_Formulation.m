classdef Dynamic_Formulation
    properties
        drive_train;
        n_DOF; % number of degrees of freedom
        M;     % Inertia matrix
        K;     % Stiffness matrix
        D;     % Damping matrix
        A;     % State matrix
        c;     % Load vector
        
        M_spc; % Inertia matrix, (SParCe version)
        K_spc; % Stiffness matrix
        D_spc; % Damping matrix
        A_spc; % State matrix
        c_spc; % Load vector
    end
    
    methods
        function obj = Dynamic_Formulation(DT)
            if(~isa(DT, 'Drivetrain') || ~strcmp(superclasses(DT), 'Drivetrain'))
                error('DT should be a Drivetrain object.');
            end
            
            obj.drive_train = DT;
            
            obj.n_DOF = obj.calculate_DOF();
            
            obj.M = obj.inertia_matrix();
            obj.K = obj.stiffness_matrix();
            obj.D = obj.damping_matrix();
            obj.A = obj.state_matrix();
            obj.c = obj.load_vector();
            
            obj.M_spc = sparse(obj.M);
            obj.K_spc = sparse(obj.K);
            obj.D_spc = sparse(obj.D);
            obj.A_spc = sparse(obj.A);
            obj.c_spc = sparse(obj.c);
            
        end
        
        function tab = disp(obj)
            
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
            [MS, EV] = eig(obj.K, obj.M);
            
            EV = diag(EV); % matrix to vector
            omega_n = sqrt(EV); % lambda to omega
            f_n = omega_n./(2.0*pi); % rad/s to Hz
            
            % sorting in ascending order:
            [f_n, idx] = sort(f_n);
            MS = MS(:, idx);
            
            % Normalizing the mode shapes so that the maximum is always +1:
            for idx = 1:length(f_n)
                [~, n] = max(abs(MS(:, idx)));

                MS(:, idx) = MS(:, idx)./MS(n, idx);
            end
            
            flag_RB = any(abs(f_n) < 1.0e-2);
            if(flag_RB)
                warning('Dynamic_Formulation:RB', 'Rigid body behavior.');
                idx = (abs(f_n) < 1.0e-2);
                f_n(idx) = 0.0;
                
                f_n = [f_n(~idx);
                       f_n(idx)];
                   
                MS = [MS(:, ~idx), MS(:, idx)];
            end
            
            %% Damped case:
            [MS_d, EV_d] = eig(obj.A);
            
            real_part = (EV_d  + EV_d')/2.0;
            imag_part = (EV_d' - EV_d )*sqrt(-1.0)/2;
            
            real_part = diag(real_part); % matrix to vector
            imag_part = diag(imag_part);
            
            omega_n = sqrt(real_part.^2 + imag_part.^2); % lambda to omega_n
            f_nd = omega_n./(2.0*pi); % rad/s to Hz
            zeta = -real_part./omega_n;
            
            % sorting in ascending order:
            [f_nd, idx] = sort(f_nd);
            zeta = zeta(idx);
            MS_d = MS_d(:, idx);
            
            % eliminating repeated values:
            f_nd = f_nd(1:2:end);
            zeta = zeta(1:2:end);
            n = obj.n_DOF(end);
            MS_d = MS_d(1:n, 1:2:end);
            
            flag_RB = any(abs(f_nd) < 1.0e-2);
            if(flag_RB)
                warning('Dynamic_Formulation:RB', 'Rigid body behavior.');
                idx = (abs(f_nd) < 1.0e-2);
                f_nd(idx) = 0.0;
                
                f_nd = [f_nd(~idx);
                        f_nd(idx)];
                   
                MS_d = [MS_d(:, ~idx), MS_d(:, idx)];
            end
            
        end
        
        function [M_r, K_r, D_r, c_r, MS_r] = modal_reduction(obj, N, prec)
            [~, MS] = obj.modal_analysis();
            [M_r, K_r, D_r, c_r, MS_r] = ...
                Dynamic_Formulation.modal_truncation(obj.M, obj.K, obj.D, obj.c, MS, N, prec);
            
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
        
        function cc = load_vector(obj)
            cc = eye(obj.n_DOF(end));
        end
        
        function nn = calculate_DOF(obj)
            nn = 2;
        end
        
        function AA = state_matrix(obj)
            n = obj.n_DOF(end);
            AA = [ zeros(n),     eye(n);
                  -obj.M\obj.K, -obj.M\obj.D];
        end
    end
    
    methods(Static)
        function [M_r, K_r, D_r, c_r, MS_r] = modal_truncation(MM, KK, DD, cc, MS, N, prec)
            MS_r = MS(:, N);
            
            eps = power(10.0, -prec);
            
            M_r = MS_r'*MM*MS_r;   M_r(abs(M_r) < eps) = 0.0;
            K_r = MS_r'*KK*MS_r;   K_r(abs(K_r) < eps) = 0.0;
            D_r = MS_r'*DD*MS_r;   D_r(abs(D_r) < eps) = 0.0;
            c_r = MS_r'*cc;        c_r(abs(c_r) < eps) = 0.0;
            
        end
                
        function time_response(A, b, T_A, T_G)
        end
        
        function [sol, time, pos, vel, acc] = Newmark(time, x0, v0, M, C, K, R, option)
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
            
            K_hat = a_1*M + a_2*C + K;
            
            %% Initial calculation:
            idx = 1;
            pos(:, idx) = x0;
            vel(:, idx) = v0;
            acc(:, idx) = M\(R(:, idx) - C*vel(:, idx) - K*pos(:, idx));
            
            %% for each time step:
            nt = nt - 1;
            
            for idx = 1:nt
                Delta_R_hat = R(:, idx + 1) + M*(a_1*pos(:, idx) + a_3*vel(:, idx) + a_4*acc(:, idx)) + ...
                                              C*(a_2*pos(:, idx) + a_5*vel(:, idx) + a_6*acc(:, idx));
                
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
        
        function [sol, time, pos, vel, acc] = Bathe(time, x0, v0, M, C, K, R)
            %BATHE solution of linear mechanical systems using Bathe's
            % method according to [1] Table 9.4.
            %
            % Parameters:
            % time: time interval for integration
            % x0: initial position
            % v0: initial velocity
            % M: Inertia matrix
            % C: Damping matrix
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
            
            K_hat1 = K + a_1*M + a_2*C;
            K_hat2 = K + a_3*M + a_4*C;
            
            R_mid = (R(:, 1:end-1) + R(:, 2:end))/2.0;
            
            %% Initial calculation:
            idx = 1;
            pos(:, idx) = x0;
            vel(:, idx) = v0;
            acc(:, idx) = M\(R(:, idx) - C*vel(:, idx) - K*pos(:, idx));
            
            %% for each time step:
            nt = nt - 1;
            
            for idx = 1:nt
                % First sub-step:
                R_hat1 = R_mid(:, idx) + M*(a_1*pos(:, idx) + a_5*vel(:, idx) + acc(:, idx)) + ...
                                         C*(a_2*pos(:, idx) +     vel(:, idx));
                pos_1 = K_hat1\R_hat1;
                vel_1 = a_2*(pos_1 - pos(:, idx)) - vel(:, idx);
%                 acc_1 = a_2*(vel_1 - vel(:, idx)) - acc(:, idx);
                
                % Second sub-step:
                R_hat2 = R(:, idx + 1) + M*(a_6*pos_1 + a_7*pos(:, idx) + a_2*vel_1 + a_8*vel(:, idx)) + ...
                                         C*(a_2*pos_1 + a_8*pos(:, idx));
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
            
            C = zeros(n);
            
            r = [0.0, 10.0]';
            
            Ts = 0.28;
            t = 0.0:Ts:(12*Ts);
            x0 = zeros(n, 1);
            v0 = x0;
            
            R = repmat(r, 1, length(t));
            
            option = 'average';
            [~, time, pos, ~, ~] = Dynamic_Formulation.Newmark(t, x0, v0, M, C, K, R, option);
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
            
            C = zeros(n);
            
            r = [0.0, 10.0]';
            
            Ts = 0.28;
            t = 0.0:Ts:(12*Ts);
            x0 = zeros(n, 1);
            v0 = x0;
            
            R = repmat(r, 1, length(t));
            
            [~, time, pos, ~, ~] = Dynamic_Formulation.Bathe(t, x0, v0, M, C, K, R);
            
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
