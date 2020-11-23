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
Gu = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_Lon, B_Lon, C(2,:), D, B_num);
Ga = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_Lon, B_Lon, C(3,:), D, B_num);
Gq = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_Lon, B_Lon, C(4,:), D, B_num);
Gt = minreal(zpk(tf(num, den)));

Gvs = V_trim * (Gt - Ga);


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
disp('-------- Q3 --------');
% inner loop TF
s = tf('s');

% q controller (inner loop)
K1 = -2*(s+0.05) / (s * (s + 0.01) ); % lead/lag
%K1 = -6.1541*(1+0.0064*s)*(1+0.035*s) / s; % PID
%K1 = -1.6667*(s+0.03)/s^2; % PI with integrator

% vs controller (outer loop)
K2 = 0.015; % proportional - works with lead/lag

G1 = Gq;
G2 = Gvs;

% inner loop sensitivity functins (de to q)
T1 = G1 * K1;
S1 = 1 / (1 + T1);
C1 = T1 / (1 + T1); 

Gin = K1 * S1 * G2;
Gin = minreal(zpk(Gin));

%sisotool(Gin)

% outer loop sensitivity functions (q to vs)
T2 = Gin * K2;
S2 = 1 / (1 + T2);
C2 = T2 / (1 + T2);


% time response
dt = 0.1;
t = 0:dt:10;

q_innerLead  = step(C1, t); 
q_inner_c = ones(1,length(t));


K_PI = -1.6667*(s+0.03)/s^2;
%K_PD = 
K_PID = -6.1541*(1+0.0064*s)*(1+0.035*s) / s;
q_innerPI = step(G1 * K_PI / (1 + G1 * K_PI), t);
q_innerPD = step(G1 * K_PI / (1 + G1 * K_PI), t);
q_innerPID = step(G1 * K_PID / (1 + G1 * K_PID), t);

figure(3); % response of different q controllers
plot(t,rad2deg(q_inner_c),'--', t,rad2deg(q_innerLead), t,rad2deg(q_innerPI), t,rad2deg(q_innerPID));
xlabel('Time (s)', 'interpreter','latex');
ylabel('Pitch Rate (deg/s)', 'interpreter','latex');
legend({'Command', 'Lead/Lag', 'PI', 'PID'}, 'interpreter','latex');
grid minor


vsc = 2.54 * ones(1,length(t)); % 500 ft/min vertical speed step input (2.54 m/s)

vs = step(C2 * 2.54, t);
qc = step(K2 * S2 * 2.54, t); 
%q  = step(C1, t); % q to qc - using step input (wrong)
%de = step(K1 * S1, t); % de to qc - using step input (wrong)
q = step(C1 * K2 * S2 * 2.54, t); % q to vs_c
de = step(K1 * S1 * K2 * S2 * 2.54, t); % de to vs_c


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
fprintf('Gust PSD\nMax ug: %.4f, std dev: %.4f\n', 3 * sqrt(meanSquare), sqrt(meanSquare)); 

meanSquare = 1/pi *trapz(omega, PSD_wg);
fprintf('Max wg: %.4f, std dev: %.4f\n\n', 3 * sqrt(meanSquare), sqrt(meanSquare)); 




%%% Gust TF
Cg = [1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 0 0 0];
Gamma_Lon = S*Gamma*Cg'; % size 5 x 3

Cq  = [0 0 1 0];
Cvs = 1/V_trim * [0 1 0 1];
D = [0 0 0];

[num, den] = ss2tf(A_Lon, Gamma_Lon, Cq, D, 2);
ug_q = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_Lon, Gamma_Lon, Cq, D, 2);
wg_q = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_Lon, Gamma_Lon, Cvs, D, 2);
ug_vs = minreal(zpk(tf(num, den)));

[num, den] = ss2tf(A_Lon, Gamma_Lon, Cvs, D, 2);
wg_vs = minreal(zpk(tf(num, den)));

%%% Max state response
% closed loop
% T = Gin*K2
% S = 1/(1+T)
% C = T/(1+T)
[mag, ~, ~] = bode(ug_q * S2, omega);
PSD_ug_q_cl = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_q_cl);
fprintf('Closed\nMax q to ug: %.4f deg/s, std dev: %.4f\n', rad2deg(3 * sqrt(meanSquare)), rad2deg(sqrt(meanSquare))); 

[mag, ~, ~] = bode(wg_q * S2, omega);
PSD_wg_q_cl = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_q_cl);
fprintf('Max q to wg: %.4f deg/s, std dev: %.4f\n', rad2deg(3 * sqrt(meanSquare)), rad2deg(sqrt(meanSquare))); 

[mag, ~, ~] = bode(ug_vs * S2, omega);
PSD_ug_vs_cl = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_vs_cl);
fprintf('Max vs to ug: %.4f m/s, std dev: %.4f\n', 3 * sqrt(meanSquare), sqrt(meanSquare)); 

[mag, ~, ~] = bode(wg_vs * S2, omega);
PSD_wg_vs_cl = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_vs_cl);
fprintf('Max vs to wg: %.4f m/s, std dev: %.4f\n\n', 3 * sqrt(meanSquare), sqrt(meanSquare)); 



% open loop 
% S = 1 so no impact
% C = Gin*K2
[mag, ~, ~] = bode(ug_q, omega);
PSD_ug_q_ol = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_q_ol);
fprintf('Open\nMax q to ug: %.4f deg/s, std dev: %.4f\n', rad2deg(3 * sqrt(meanSquare)), rad2deg(sqrt(meanSquare))); 

[mag, ~, ~] = bode(wg_q, omega);
PSD_wg_q_ol = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_q_ol);
fprintf('Max q to wg: %.4f deg/s, std dev: %.4f\n', rad2deg(3 * sqrt(meanSquare)), rad2deg(sqrt(meanSquare))); 

[mag, ~, ~] = bode(ug_vs, omega);
PSD_ug_vs_ol = (squeeze(mag.^2))' .* PSD_ug;
meanSquare = 1/pi *trapz(omega, PSD_ug_vs_ol);
fprintf('Max vs to ug: %.4f m/s, std dev: %.4f\n', 3 * sqrt(meanSquare), sqrt(meanSquare)); 

[mag, ~, ~] = bode(wg_vs, omega);
PSD_wg_vs_ol = (squeeze(mag.^2))' .* PSD_wg;
meanSquare = 1/pi *trapz(omega, PSD_wg_vs_ol);
fprintf('Max vs to wg: %.4f m/s, std dev: %.4f\n', 3 * sqrt(meanSquare), sqrt(meanSquare)); 


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



G_ug_de = 1/Gq * 1/ug_q;
G_wg_de = 1/Gq * 1/wg_q;
G_ug_de1 = 1/Gvs * 1/ug_vs;
G_wg_de1 = 1/Gvs * 1/wg_vs;







%% Q4



