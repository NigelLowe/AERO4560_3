%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main Aircraft Simulation Script File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD FLIGHT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ FlightData ] = LoadFlightData_aircraft4_90kts_CG1();
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load X0 and U0 Here!
load ICs_aircraft4_90Kts_CG1.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION VECTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DT = 0.01;          % Integration time interval
T0 = DT;            % Simulation start time
TF = 10;            % Termination time for simulation
    
n_pts = round((TF-T0)/DT+1);
X = zeros(12,n_pts);
U = zeros(5,n_pts);
T = zeros(1,n_pts);

X(:,1) = X0;
U(:,1) = U0;
T(1)   = DT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIM CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_trim  = sqrt(X0(1)^2+X0(2)^2+X0(3)^2);     % (m/s)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:n_pts % Start Simulation loop
    
    % Current Time
    T(i) = i*DT;

	% Control Input
    U(:,i-1)=U0;
    
    if T(i)>=1
        U(:,i-1)=U0+[0;1*pi/180;0;0;0];
    end
     
    % Gust Input ([u,v,w,p,q,r]^T gust components)
    Xg = [0;0;0;0;0;0];
    
    % Integration
    [X_out] = aero4560_euler(DT,X(:,i-1),Xg,U(:,i-1), FlightData);
    X(:,i) = X_out;

end     % End Simulation loop

U(:,i) = U(:,i-1); % Save final control inputs

figure;plot(T,X(1,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SIMULATION RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
