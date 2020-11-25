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


cd(Current_Folder); % go back inside longitudinal folder

% longitudinal selection matrix
G = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0];
A_long = G*A*G'; % u, w, q, theta

Akin_Lon = G*Akin*G';

% longitudinal control selection
H = [1 0 0 0 0;
     0 1 0 0 0];
B_long = G*B*H'; % dT, de

C = [1 0 0 0;
     0 1/V_trim 0 0;
     0 0 1 0;
     0 0 0 1];
D = [0 0];  


[num, den] = ss2tf(A_long, B_long, C(1,:), D, 1);
G_dt_u = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, B_long, C(2,:), D, 1);
G_dt_a = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, B_long, C(3,:), D, 1);
G_dt_q = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, B_long, C(4,:), D, 1);
G_dt_t = minreal(zpk(tf(num, den)));

G_dt_vs = V_trim * (G_dt_t - G_dt_a);



[num, den] = ss2tf(A_long, B_long, C(1,:), D, 2);
G_de_u = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, B_long, C(2,:), D, 2);
G_de_a = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, B_long, C(3,:), D, 2);
G_de_q = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, B_long, C(4,:), D, 2);
G_de_t = minreal(zpk(tf(num, den)));

G_de_vs = V_trim * (G_de_t - G_de_a);


%sisotool(Gq);

%{
%% Q1
Use root locus tools to investigate why the design of an autopilot based on
using elevator to directly control vertical speed is problematic. Discuss
the nature of the inherent problems. (500ft/min vertical speed step input.
Rise time 3s, settling time 10s)
500ft/min = 2.54 m/s

Zeros are one large stable, one large unstable, one small stable root.
Increasing gain to control vertical speed will push the small stable root
into the right hand plane 



Proportional controller - improve rise time
Derivative controller - reduce overshoot
Integral controller - reduce steady state error



       rise time    overshoot    settling time    steady state error
Kp         v            ^            ~                     v               
Ki         v            ^            ^                     v
Kd         ~            v            v                     -
lead/lag   v            v                                  v
v = decrease, ^ = increase, ~ = small change, - = no change


C = -0.0001, pole at 30.04: Trying to cancel the unstable zero is very difficult becuase it has to be
exactly the same value root to cancel but they are defined by a lot of
decimal places. So the unstable zero will exist in the root locus even when
you try cancel it. 


Proportional C = 0.0001: stable response. Rise time = 3.62s, peak amplitude
= -910 (58.5% overshoot) at 11.3s, settling time = 54.1s, steady state = -573 
Decrease to C = 1e-10: tr = 3.42s, peak -853 (57.3%) at 10.9s, ts = 51.8,
ss = -542. Proportional controller by itself is not enough to create the
desired response.

Other way of decreasing rise time by adding integrator but increases
overshoot and settling time. 








Q2
*PI with integrator - within command bandwidth(integrator to get rid of differentiator in tf)
*PI - only way to get low sse is with extremely fast response. < 1
millisecond but cross over frequency at 1000 rad/s (has 10 times gain as PI
with integrator)
*PI - got cross over frequency 2-10 rad/s with -0.2875(s+0.032)/s. But sse
is no where near desired (at 0.1 amplitude)
*Lead/Lag - want pole to be as close to zero as possible to be like a
second integrator for the sse response but needs be be within 100 times of zero.

converts commanded vertical speed to a pitch rate

%}

%% Q3
close all

disp('-------- Q3 --------');
% inner loop TF
s = tf('s');

% q controller (inner loop)
K1_q = -2*(s+0.5) / (s * (s + 0.01) ); % lead/lag - cross over 9.4 rad/s
%K1_q = -6.1541*(1+0.0064*s)*(1+0.035*s) / s; % PID
%K1_q = -1.6667*(s+0.03)/s^2; % PI with integrator

% vs controller (outer loop)
%K2_vs = 0.015; % proportional 
K2_vs = 0.03 / (s+3); % with pole - least sudden response
%K2_vs = 0.027404*(s+0.1918)/s; % PI
G1_vs = G_de_q;
G2_vs = G_de_vs;

% inner loop sensitivity functins (de to q)
T1_vs = G1_vs * K1_q;
S1_vs = 1 / (1 + T1_vs);
C1_vs = T1_vs / (1 + T1_vs); 

Gin_vs = K1_q * S1_vs * G2_vs;
Gin_vs = minreal(zpk(Gin_vs));

%sisotool(Gin)

% outer loop sensitivity functions (q to vs)
T2_vs = Gin_vs * K2_vs;
S2_vs = 1 / (1 + T2_vs);
C2_vs = T2_vs / (1 + T2_vs);


% time response
dt = 0.01;
t = 0:dt:10;

%q_innerLead  = step(C1, t); 
q_innerLead  = step(C1_vs, t); 
q_inner_c = ones(1,length(t));


% other controllers for comparison
K_P = -27; % cross over 1500 rad/s
K_PI = -0.0073788*(s+285.4)/s; % cross over 9.7 rad/s
K_PI_mod = -1.6667*(s+0.03)/s^2; % cross over 8.2 rad/s
K_PD = -52.723;
K_PID = -6.1541*(1+0.0064*s)*(1+0.035*s) / s; % cross over 38.5 rad/s
q_innerPI = step(G1_vs * K_PI / (1 + G1_vs * K_PI), t);
q_innerPI_mod = step(G1_vs * K_PI_mod / (1 + G1_vs * K_PI_mod), t);
q_innerPD = step(G1_vs * K_PD / (1 + G1_vs * K_PD), t);
q_innerPID = step(G1_vs * K_PID / (1 + G1_vs * K_PID), t);

plotXmax = 5;
plotYmin = 40;
plotYmax = 65;
figure(3); % response of different q controllers
plot(t,rad2deg(q_inner_c),'--', t,rad2deg(q_innerLead), t,rad2deg(q_innerPI), t,rad2deg(q_innerPI_mod), t,rad2deg(q_innerPID)); 
hold on
patch([3 plotXmax plotXmax 3], [plotYmin plotYmin 0.98 * rad2deg(1) 0.98 * rad2deg(1)], [0.8 0.8 0.8])
patch([3 plotXmax plotXmax 3], [1.02 * rad2deg(1) 1.02 * rad2deg(1) plotYmax plotYmax], [0.8 0.8 0.8])
xlabel('Time (s)', 'interpreter','latex');
ylabel('Pitch Rate (deg/s)', 'interpreter','latex');
legend({'Command', 'Lag', 'PI', 'PI Modified', 'PID'}, 'interpreter','latex');
ylim([plotYmin plotYmax])
xlim([0 5])
grid minor


% step input to vertical speed
vsc = 2.54 * ones(1,length(t)); % 500 ft/min vertical speed step input (2.54 m/s)

vs = step(C2_vs * 2.54, t);
qc = step(K2_vs * S2_vs * 2.54, t); 
q = step(C1_vs * K2_vs * S2_vs * 2.54, t); % q to vs_c
de = step(K1_q * S1_vs * K2_vs * S2_vs * 2.54, t); % de to vs_c


figure(4);
plot(t,vs, t,vsc);
xlabel('Time (s)', 'interpreter','latex');
ylabel('Vertical Speed (m/s)', 'interpreter','latex');
legend({'$v_s$', '$v_{s,c}$'}, 'interpreter','latex');
grid minor

figure(5);
plot(t,rad2deg(q), t,rad2deg(qc));
xlabel('Time (s)', 'interpreter','latex');
ylabel('Pitch Rate (deg/s)', 'interpreter','latex');
legend({'$q$', '$q_c$'}, 'interpreter','latex');
grid minor

trim_de = U0(2);
figure(6);
plot(t,rad2deg(de+trim_de));
xlabel('Time (s)', 'interpreter','latex');
ylabel('Elevator deflection (deg)', 'interpreter','latex');
grid minor

fprintf('From plot\nMax vs: %.4f m/s\n', max(vs)); 
fprintf('Max q: %.4f deg/s\n', max(rad2deg(q))); 
fprintf('Max de: %.4f, min de: %.4f deg\n\n', max(rad2deg(de+trim_de)), min(rad2deg(de+trim_de))); 



% Gust input

% Ref wind at 20 feet AGL
u20 = convvel(35,'kts','m/s'); 
h_ref = convlength(250, 'ft','m'); 

sigma_w = 0.1*u20;
sigma_u = sigma_w/(0.177+0.0027*h_ref)^0.4;
sigma_v = sigma_u;

% Set length scales from slides
Lw = h_ref;
Lu = h_ref/(0.177+0.0027*h_ref)^1.2;
Lv = Lu;

b = FlightData.Geo.b;
omega = logspace(-5, 5, 10000); 

%%% PSD
PSD_ug = sigma_u.^2.*2.*Lu./pi./V_trim./(1+(Lu./V_trim.*omega).^2);
PSD_wg = sigma_w.^2*Lw./pi./V_trim.*(1+3.*(Lw./V_trim.*omega).^2)./((1+(Lw./V_trim.*omega).^2).^2);
PSD_qg = (omega./V_trim).^2./(1+(4.*b./pi./V_trim.*omega).^2).*PSD_wg;

figure(7);
subplot(1,3,1)
loglog(omega, PSD_ug);
xlabel('Frequency \omega (rad/s)');
ylabel('$\Phi_{u_g}$', 'interpreter','latex');

subplot(1,3,2)
loglog(omega, PSD_wg);
xlabel('Frequency \omega (rad/s)');
ylabel('$\Phi_{w_g}$', 'interpreter','latex');

subplot(1,3,3)
loglog(omega, PSD_qg);
xlabel('Frequency \omega (rad/s)');
ylabel('$\Phi_{q_g}$', 'interpreter','latex');

%%% Max gust
meanSquare = 1/pi *trapz(omega, PSD_ug);
fprintf('Gust PSD\nMax ug: %.4f\n', 3 * sqrt(meanSquare)); 

meanSquare = 1/pi *trapz(omega, PSD_wg);
fprintf('Max wg: %.4f\n\n', 3 * sqrt(meanSquare)); 




%%% Gust TF
Cg = [1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 0 0 0];
Gamma_Lon = G*Gamma*Cg'; % size 5 x 3

Cq  = [0 0  1 0];
Cvs = [0 -1 0 V_trim];
D = [0 0 0];

[num, den] = ss2tf(A_long, Gamma_Lon, Cq, D, 1);
q_ug = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, Gamma_Lon, Cq, D, 2);
q_wg = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, Gamma_Lon, Cvs, D, 1);
vs_ug = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_long, Gamma_Lon, Cvs, D, 2);
vs_wg = minreal(zpk(tf(num, den)));

%%% Max state response
% closed loop
% T = Gin*K2
% S = 1/(1+T)
% C = T/(1+T)
[mag, ~, ~] = bode(q_ug * S2_vs, omega);
PSD_ug_q_cl = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_q_cl);
fprintf('Closed\nMax q to ug: %.4f deg/s\n', rad2deg(3 * sqrt(meanSquare))); 

[mag, ~, ~] = bode(q_wg * S2_vs, omega);
PSD_wg_q_cl = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_q_cl);
fprintf('Max q to wg: %.4f deg/s\n', rad2deg(3 * sqrt(meanSquare))); 

[mag, ~, ~] = bode(vs_ug * S2_vs, omega);
PSD_ug_vs_cl = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_vs_cl);
fprintf('Max vs to ug: %.4f m/s\n', 3 * sqrt(meanSquare)); 

[mag, ~, ~] = bode(vs_wg * S2_vs, omega);
PSD_wg_vs_cl = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_vs_cl);
fprintf('Max vs to wg: %.4f m/s\n\n', 3 * sqrt(meanSquare)); 



% open loop 
% S = 1 so no impact
% C = Gin*K2
[mag, ~, ~] = bode(q_ug, omega);
PSD_ug_q_ol = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_q_ol);
fprintf('Open\nMax q to ug: %.4f deg/s\n', rad2deg(3 * sqrt(meanSquare))); 

[mag, ~, ~] = bode(q_wg, omega);
PSD_wg_q_ol = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_q_ol);
fprintf('Max q to wg: %.4f deg/s\n', rad2deg(3 * sqrt(meanSquare))); 

[mag, ~, ~] = bode(vs_ug, omega);
PSD_ug_vs_ol = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_vs_ol);
fprintf('Max vs to ug: %.4f m/s\n', 3 * sqrt(meanSquare)); 

[mag, ~, ~] = bode(vs_wg, omega);
PSD_wg_vs_ol = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_vs_ol);
fprintf('Max vs to wg: %.4f m/s\n\n', 3 * sqrt(meanSquare)); 


figure(8);
subplot(2,2,1)
loglog(omega, PSD_ug_q_cl);
hold on
loglog(omega, PSD_ug_q_ol,'--');
xlabel('Frequency \omega (rad/s)');
ylabel('$\Phi_{u_g, q}$', 'interpreter','latex');
%ylim([10^(-40) 10^2])

subplot(2,2,2)
loglog(omega, PSD_wg_q_cl);
hold on
loglog(omega, PSD_wg_q_ol,'--');
xlabel('Frequency \omega (rad/s)');
ylabel('$\Phi_{w_g, q}$', 'interpreter','latex');
%ylim([10^(-40) 10^2])

subplot(2,2,3)
loglog(omega, PSD_ug_vs_cl);
hold on
loglog(omega, PSD_ug_vs_ol,'--');
xlabel('Frequency \omega (rad/s)');
ylabel('$\Phi_{u_g, v_s}$', 'interpreter','latex');
%ylim([10^(-40) 10^2])

subplot(2,2,4)
loglog(omega, PSD_wg_vs_cl);
hold on
loglog(omega, PSD_wg_vs_ol,'--');
xlabel('Frequency \omega (rad/s)');
ylabel('$\Phi_{w_g, v_s}$', 'interpreter','latex');
%ylim([10^(-40) 10^2])

legend('Closed Loop', 'Open Loop');



% including each gust type - but G_gust_control dont match up between q and
% vs
G_de_ug = 1/G_de_q * q_ug * S2_vs;
G_de_wg = 1/G_de_q * q_wg * S2_vs;
G_de_ug1 = 1/G_de_vs * vs_ug * S2_vs;
G_de_wg1 = 1/G_de_vs * vs_wg * S2_vs;

%{
ug_de = step(G_ug_de, t);
ug_de1 = step(G_ug_de1, t);
wg_de = step(G_wg_de, t);
wg_de1 = step(G_wg_de1, t);

figure(9)
subplot(1,2,1)
plot(t,ug_de,'k', t,ug_de1,'--');
xlabel('t (s)');
ylabel('$u_g \delta_e$', 'interpreter','latex');
legend('q', 'vs');

subplot(1,2,2)
plot(t,wg_de,'b', t,wg_de1,'r');
xlabel('t (s)');
ylabel('$w_g \delta_e$', 'interpreter','latex');
%}
[mag, ~, ~] = bode(G_de_ug, omega);
PSD_de_ug = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_de_ug);
fprintf('Max de (w/ q) to ug: %.4f deg\n', rad2deg(3 * sqrt(meanSquare))); 

[mag, ~, ~] = bode(G_de_wg, omega);
PSD_de_wg = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_de_wg);
fprintf('Max de (w/ q) to wg: %.4f deg\n', rad2deg(3 * sqrt(meanSquare))); 

[mag, ~, ~] = bode(G_de_ug1, omega);
PSD_de_ug1 = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_de_ug1);
fprintf('Max de (w/ vs) to ug: %.4f deg\n', rad2deg(3 * sqrt(meanSquare))); 

[mag, ~, ~] = bode(G_de_wg1, omega);
PSD_de_wg1 = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_de_wg1);
fprintf('Max de (w/ vs) to wg: %.4f deg\n', rad2deg(3 * sqrt(meanSquare))); 


%{
G_q_d1 = (1+G1*K1+G2*K1*K2+G1*G2*K1^2*K2) / ((1+G1*K1)*(1+G1*K1+G2*K1*K2));
G_q_d2 = -G1*K1*K2 / (1+G1*K1+G2*K1*K2);
G_vs_d1 = -G2*K1 / (1+G1*K1+G2*K1*K2);
G_vs_d2 = (1+G1*K1) / (1+G1*K1+G2*K1*K2);
%}

figure(9);
subplot(1,2,1)
loglog(omega, PSD_de_ug);
hold on
loglog(omega, PSD_de_ug1,'--');
xlabel('Frequency \omega (rad/s)');
ylabel('$u_g$', 'interpreter','latex');
%ylim([10^(-40) 10^2])

subplot(1,2,2)
loglog(omega, PSD_de_wg);
hold on
loglog(omega, PSD_de_wg1,'--');
xlabel('Frequency \omega (rad/s)');
ylabel('$w_g$', 'interpreter','latex');
%ylim([10^(-40) 10^2])

legend('with q', 'with v_s');



%% Q4
%{

% Make the closed loop K
% K_PID = tf([12 12.012 0.012],[1 0]); % Test PID
K_PID = tf([-0.12, -0.8, -8],[1 0]);
K_PID_3 = tf([-0.01, -384,-384],[1,0]);
    
% K_PID = K_PID_3;
sys_cl = K_PID*G_q_de/(1+K_PID*G_q_de);

% asd = zpk([-0.5,-0.52],[0],-1);

% sys_cl_2 = asd*G_q_de/(1+asd*G_q_de);

% figure()
% step(G_q_de,0:0.05:20)
% figure()
% step(sys_cl,0:0.05:20)
% % figure()
% % step(sys_cl_2,0:0.05:20)
% figure()
% bode(G_q_de*K_PID)
% margin(K_PID*G_q_de)


% G_dt_du
 G_dt_du = -G_dt_vs * K2 * S_vs * K_q * S_q * G_de_u ;

% Plant Transfer Function for u - loop
 G_dtvs_u = G_dt_u + G_dt_du ;
 K_u = 1;
 S_u = 1 / (1 + K_u * G_dt_u );


 % Re design Vs loop
 G_de_dvs = -G_de_u * K_u * S_u * G_dt_vs ;
 G_deu_vs = G_de_vs + G_de_dvs ;
 K_vs_2 = -0.0024655; 
 % I don't know if 
 % New G_dtvs_u_2
 % G_dt_du
 G_dt_du = -G_dt_vs * K_vs_2 * S_vs * K_q * S_q * G_de_u ;
 % Plant Transfer Function for u - loop
 G_dtvs_u_2 = G_dt_u + G_dt_du ;
 
 
 
 %}
K_u = 1;

%% Q5
DT = 0.05;           % Integration time interval
T0 = DT;            % Simulation start time
TF = 15;           % Termination time for simulation
    
n_pts = round((TF-T0)/DT+1);
T = zeros(1,n_pts);
T(1)   = DT;

clear X U
X(:,1) = G * zeros(12,1); %X0;
U(:,1) = H * zeros(5,1); %U0;

Y(:,1) = zeros(5,1);
D = 0;

Cvs = [1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1;
       0 -1 0 V_trim];
   
Dvs = [0 0];
[Ao_vs, Bo_vs, Co_vs, Do_vs] = ssdata(ss(K2_vs)); % outer loop - vertical speed
[Ai_q, Bi_q, Ci_q, Di_q] = ssdata(ss(K1_q)); % inner loop - pitch rate
[Au , Bu , Cu , Du] = ssdata (ss( K_u ));

Xo = zeros(size(Ao_vs,1),n_pts);
Xi = zeros(size(Ai_q,1),n_pts); 

vsc = 2.54 * ones(1,n_pts);
for i = 2:n_pts % Start Simulation loop

    % Current Time
    T(i) = i*DT;
    
    % Gust Input ([u,v,w,p,q,r]^T gust components)
    Xg = [0;0;0;0;0;0];
    
    
    if T(i) > 1
        % add control loops
        % vertical speed guidance loops
        % vs - outer loop
        eo = vsc(i-1) - Y(5,i-1);
        Xo(:,i) = aero4560_euler(DT, Ao_vs, Bo_vs, Xo(:,i-1), eo);
        ic = Co_vs * Xo(:,i) + Do_vs * eo; % [vs; qc]

        % q - inner loop
        ei = ic - Y(3,i-1); % Xo(:,i-1); % Xo(1,i) = q
        Xi(:,i) = aero4560_euler(DT, Ai_q, Bi_q, Xi(:,i-1), ei);
        de = Ci_q * Xi(:,i) + Di_q * ei;


        % auto-throttle %%% add you stuff here
        % u - air speed
        % dT = 

        U(:,i) = [0; de];
        %U(:,i) = [dT; de]; % after you add q4 loop
    else
         U(:,i) = [0; 0];
    end
    
    Xdot = A_long * X(:,i-1) + B_long * U(:,i);
    X(:,i) = X(:,i-1) + Xdot*DT;
   
    Y(:,i) = Cvs * X(:,i) + Dvs * U(:,i); % [u, alpha, q, theta, vs]   
    
end     % End Simulation loop


figure(10); % linear
subplot(1,3,1)
plot(T, Y(5,:), T,vsc,'--');
ylabel('$v_s$ (m/s)', 'interpreter','latex');
xlabel('Time (s)');
legend({'$v_s$', '$v_{s,c}$'}, 'interpreter','latex');
xlim([0 10]);
grid minor

subplot(1,3,2)
plot(T, rad2deg(Y(3,:)+X0(5)));
ylabel('$q$ (deg/s)', 'interpreter','latex');
xlabel('Time (s)');
xlim([0 10]);
grid minor

subplot(1,3,3)
plot(T, rad2deg(U(2,:)+U0(2)));
ylabel('$\delta_e$ (deg)', 'interpreter','latex');
xlabel('Time (s)');
xlim([0 10]);
grid minor

