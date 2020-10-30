% Assignment 3 - AERO4560
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Main Aircraft Simulation Script File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

set(groot,'defaultLineLineWidth',2.0,...
    'DefaultAxesFontSize', 20, ...
    'defaultLineMarkerSize',30,...
    'defaultAxesXGrid','on',...
    'defaultAxesYGrid','on')

loadCase = 1; 
% 1 = 50kts, CG1
% 2 = 50kts, CG2
% 3 = 90kts, CG1
% 4 = 90kts, CG2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD FLIGHT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Current_Folder = pwd;
TopFolder = fileparts(pwd);
cd(TopFolder);

switch loadCase
    case 1
        [ FlightData ] = LoadFlightData_aircraft4_90kts_CG1();
        load ICs_aircraft4_90Kts_CG1.mat

    case 2
        [ FlightData ] = LoadFlightData_aircraft4_90kts_CG2();
        load ICs_aircraft4_90Kts_CG2.mat
        
    case 3
        [ FlightData ] = LoadFlightData_aircraft4_50kts_CG1();
        load ICs_aircraft4_50Kts_CG1.mat

    case 4
        [ FlightData ] = LoadFlightData_aircraft4_50kts_CG2();
        load ICs_aircraft4_50Kts_CG2.mat
        
end

%cd(Current_Folder); % go back inside longitudinal folder


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TRIM CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_trim  = sqrt(X0(1)^2+X0(2)^2+X0(3)^2);     % (m/s)
alt0    = -X0(12); % (m)
alt0ft  = alt0/0.3048;        % (ft)

%{
disp(' Trimmed State and Controls')
long_trim = ([' V = ',num2str(V_trim),' m/s, ',' Alpha = ',num2str(X0(3)/V_trim*57.3),' deg, ',' Theta = ',num2str(X0(8)*57.3),' deg. ']);
lat_trim  = ([' Beta = ',num2str(X0(2)/V_trim*57.3),' deg. ',' Phi = ',num2str(X0(7)*57.3),' deg, ',' Psi = ',num2str(X0(9)*57.3),' deg. ']);
Cont_trim = ([' Throttle = ',num2str(U0(1)),' (0-1) ',' Elevator = ',num2str(U0(2)*57.3),' deg, ',' Aileron = ',num2str(U0(3)*57.3),' deg, ',' Rudder = ',num2str(U0(4)*57.3),' deg. ',' Flap = ',num2str(U0(5)*57.3),' deg. ']);
disp(long_trim)
disp(lat_trim)
disp(Cont_trim)
%}
        
 %%%%%%%%%%%%%%%%%%%%%%%%% DETERMINE SYSTEM MATRICES %%%%%%%%%%%%%%%
% Iteration:
Xdot = zeros(12,1);
delta = 1;
while delta > 1e-8
    XdotPrev = Xdot;
    [ForceCoeff0, MomentCoeff0] = aero4560_aero(X0, zeros(6,1), Xdot, U0, FlightData);
    [Xdot] = aero4560_motion(X0, ForceCoeff0, MomentCoeff0, FlightData);
    delta = sumabs(Xdot - XdotPrev);
end
Xdot0 = Xdot;

% Create A matrix by perturbing states
perturb = 1e-7;
for i = 1:12
    Xa = X0;
    Xa(i) = X0(i) + perturb;
    Xdot = zeros(12,1);
    delta = 1;
    while delta > 1e-8
        XdotPrev = Xdot;
        [ForceCoeff, MomentCoeff] = aero4560_aero(Xa, zeros(6,1), Xdot, U0, FlightData);
        [Xdot] = aero4560_motion(Xa, ForceCoeff, MomentCoeff, FlightData);
        delta = sumabs(Xdot - XdotPrev);
    end
    A(:,i) = (Xdot - Xdot0)/perturb;
end

% Create B matrix by perturbing control inputs
for i = 1:5
    Ub = U0;
    Ub(i) = U0(i) + perturb;
    Xdot = zeros(12,1);
    delta = 1;
    while delta > 1e-8
        XdotPrev = Xdot;
        [ForceCoeff, MomentCoeff] = aero4560_aero(X0, zeros(6,1), Xdot, Ub, FlightData);
        [Xdot] = aero4560_motion(X0, ForceCoeff, MomentCoeff, FlightData);
        delta = sumabs(Xdot - XdotPrev);
    end
    B(:,i) = (Xdot - Xdot0)/perturb;
end

% Create Akin matrix by perturbing states and not updating aerodynamics

for i = 1:12
    X = X0;
    X(i) = X0(i) + perturb;
    [Xdot] = aero4560_motion(X, ForceCoeff0, MomentCoeff0, FlightData);
    
    Akin(:,i) = (Xdot - Xdot0)/perturb;
end
Gamma = Akin - A;



% longitudinal selection matrix
S = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0];
A_Lon = S*A*S'; % u, w, q, theta

Akin_Lon = S*Akin*S';

% longitudinal control selection
G = [1 0 0 0 0;
     0 1 0 0 0];
B_Lon = S*B*G'; % dT, de

C = [1 0 0 0;
     0 1/V_trim 0 0;
     0 0 1 0;
     0 0 0 1];
D = [0 0];  


B_num = 2; %% 1 fr throttle, 2 for elevator


[num, den] = ss2tf(A_Lon, B_Lon, C(1,:), D, B_num);
Gu = tf(num, den);
disp('u:');
minreal(zpk(Gu))

[num, den] = ss2tf(A_Lon, B_Lon, C(2,:), D, B_num);
Ga = tf(num, den);
disp('alpha:');
minreal(zpk(Ga))

[num, den] = ss2tf(A_Lon, B_Lon, C(3,:), D, B_num);
Gp = tf(num, den);
disp('p:');
minreal(zpk(Gp))

[num, den] = ss2tf(A_Lon, B_Lon, C(4,:), D, B_num);
Gt = tf(num, den);
disp('theta:');
minreal(zpk(Gt))

Gvs = V_trim * (Gt - Ga);
disp('vs:');
minreal(zpk(Gvs))

