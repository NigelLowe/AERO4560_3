%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main Aircraft Simulation Script File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD FLIGHT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
clear all
clc
for cases = 1:4
    switch cases
        case 1
            % Flight data for CG 1 (nominal)
            [ FlightData ] = LoadFlightData_aircraft4_90kts_CG1;
            % Initial Condition
            load ICs_aircraft4_90Kts_CG1.mat
            
        case 2
            % Flight data for CG 2
            [ FlightData ] = LoadFlightData_aircraft4_90kts_CG2;
            % Initial Condition
            load ICs_aircraft4_90Kts_CG2.mat
            
        case 3
            % Flight data for CG 1 (nominal)
            [ FlightData ] = LoadFlightData_aircraft4_50kts_CG1;
            % Initial Condition
            load ICs_aircraft4_50Kts_CG1.mat
            
        case 4
            % Flight data for CG 2
            [ FlightData ] = LoadFlightData_aircraft4_50kts_CG2;
            % Initial Condition
            load ICs_aircraft4_50Kts_CG2.mat
    end
    %%
    % Xdot initially
    % Obtain state derivatives (Iterate more)
    error = 1e-8;
    diff = 1;
    Xdot = zeros(12,1);
    while diff > error
        Xdot0 = Xdot;
        [ForceCoeff, MomentCoeff] = aero4560_aero(X0,zeros(6,1),Xdot,U0,FlightData);
        [Xdot] = aero4560_motion(X0, ForceCoeff, MomentCoeff, FlightData);
        diff = sumabs(Xdot - Xdot0);
    end
    Xdot0 = Xdot;
    %% Perturbing states to get A matrix
    dx = 1e-8;
    
    for i = 1:12
        X = X0;
        X(i) = X0(i) + dx;
        Xdot = zeros(12,1);
        diff2 = 1;
        while diff2 > error
            XdotP = Xdot;
            [ForceCoeff, MomentCoeff] = aero4560_aero(X,zeros(6,1),Xdot,U0,FlightData);
            [Xdot] = aero4560_motion(X, ForceCoeff, MomentCoeff, FlightData);
            diff2 = sumabs(Xdot - XdotP);
        end
        XdotP = Xdot;
        A(:,i) = (XdotP - Xdot0)/dx;
    end
    %% Get original X_dot
    error = 1e-8;
    diff = 1;
    Xdot = zeros(12,1);
    while diff > error
        Xdot0 = Xdot;
        [ForceCoeff, MomentCoeff] = aero4560_aero(X0,zeros(6,1),Xdot,U0,FlightData);
        [Xdot] = aero4560_motion(X0, ForceCoeff, MomentCoeff, FlightData);
        diff = sumabs(Xdot - Xdot0);
    end
    %% Perturbing states to get Akin matrix
    dx = 1e-8;
    
    for i = 1:12
        X = X0;
        X(i) = X0(i) + dx;
        Xdot = zeros(12,1);
        [Xdot] = aero4560_motion(X, ForceCoeff, MomentCoeff, FlightData);
        
        Akin(:,i) = (Xdot - Xdot0)/dx;
    end
    
    Gamma = Akin - A;
    C = [1 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0 0];
    Gamma_Lon = C*Gamma*transpose(C);
    
    C2 = [0 1 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 0 0 0];
    
    Gamma_Lat = C2*Gamma*transpose(C2);
    
    %%
    error = 1e-8;
    diff = 1;
    Xdot = zeros(12,1);
    while diff > error
        Xdot0 = Xdot;
        [ForceCoeff, MomentCoeff] = aero4560_aero(X0,zeros(6,1),Xdot,U0,FlightData);
        [Xdot] = aero4560_motion(X0, ForceCoeff, MomentCoeff, FlightData);
        diff = sumabs(Xdot - Xdot0);
    end
    Xdot0 = Xdot;
    %% B matrix
    du = 1e-8;
    for i = 1:5
        U = U0;
        U(i) = U0(i) + du;
        Xdot = zeros(12,1);
        diff3 = 1;
        while diff3 > error
            XdotP = Xdot;
            [ForceCoeff, MomentCoeff] = aero4560_aero(X0,zeros(6,1),Xdot,U,FlightData);
            [Xdot] = aero4560_motion(X0, ForceCoeff, MomentCoeff, FlightData);
            diff3 = sumabs(Xdot - XdotP);
        end
        XdotP = Xdot;
        B(:,i) = (XdotP - Xdot0)/dx;
    end
    
    
    %% Lateral matrix
    C2 = [0 1 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 0 0 0];
    A_Lat = C2*A*transpose(C2);
    
    G2 = [ 0 0;
        0 0;
        1 0;
        0 1;
        0 0;];
    B_Lat = C2*B*G2;
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIM CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    V_trim  = sqrt(X0(1)^2+X0(2)^2+X0(3)^2);     % (m/s)
    V = V_trim;
    alt0    = -X0(12); % (m)
    alt0ft  = alt0/0.3048;        % (ft)
    alpha_0 = X0(3)/V_trim; % (rad)
    beta_0  = 0/57.3; % (rad)
    gamma   = 0/57.3; % (rad)
    phi_0   = 0/57.3; % (rad)
    dp      = U0(1);     % (fraction of travel)
    de      = U0(2); % (rad)
    da      = 0/57.3; % (rad)
    dr      = 0/57.3; % (rad)
    df      = U0(5); % (rad)
    
    disp(' ')
    disp(' Trimmed State and Controls')
    long_trim = ([' V = ',num2str(V_trim),' m/s, ',' Alpha = ',num2str(X0(3)/V_trim*57.3),' deg, ',' Theta = ',num2str(X0(8)*57.3),' deg. ']);
    lat_trim  = ([' Beta = ',num2str(X0(2)/V_trim*57.3),' deg. ',' Phi = ',num2str(X0(7)*57.3),' deg, ',' Psi = ',num2str(X0(9)*57.3),' deg. ']);
    Cont_trim = ([' Throttle = ',num2str(U0(1)),' (0-1) ',' Elevator = ',num2str(U0(2)*57.3),' deg, ',' Aileron = ',num2str(U0(3)*57.3),' deg, ',' Rudder = ',num2str(U0(4)*57.3),' deg. ',' Flap = ',num2str(U0(5)*57.3),' deg. ']);
    disp(' ')
    disp(long_trim)
    disp(' ')
    disp(lat_trim)
    disp(' ')
    disp(Cont_trim)
    disp(' ')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% DETERMINE SYSTEM MATRICES %%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    C_b = [1/V 0 0 0 0];
    C_p = [0 1 0 0 0];
    C_r = [0 0 1 0 0];
    C_ph = [0 0 0 1 0];
    C_ps = [0 0 0 0 1];
    
    
    %% Transfer functions
    % Choose appropriate C for each state
    D = [0 0];
    % beta
    [n_b, d_b] = ss2tf(A_Lat, B_Lat, C_b,D,2);
    tfb = tf(n_b, d_b);
    tfb = zpk(tfb);
    % p
    [n_p, d_p] = ss2tf(A_Lat, B_Lat, C_p,D,2);
    tfp = tf(n_p, d_p);
    tfp = zpk(tfp);
    % r
    [n_r, d_r] = ss2tf(A_Lat, B_Lat, C_r,D,2);
    tfr = tf(n_r, d_r);
    tfr = zpk(tfr);
    % phi
    [n_ph, d_ph] = ss2tf(A_Lat, B_Lat, C_ph,D,2);
    tfph = tf(n_ph, d_ph);
    tfph = zpk(tfph);
    % psi
    [n_ps, d_ps] = ss2tf(A_Lat, B_Lat, C_ps,D,2);
    tfps = tf(n_ps, d_ps);
    tfps = zpk(tfps);
end


