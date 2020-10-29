%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates and returns the flight data structure which
% defines the aircraft data and aerodynamic derivatives.
%
% Aircraft data and aerodynamic derivatives  - CG 22% mac.
% Flight config = clean_xx - flaps deflected 0 deg, clean configuration
% 
% (c) Peter W. Gibbens & S. Dumble, 1 March, 2011.
%
% Updates 8 May, 2014 to include flap derivatives
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ FD ] = Spero_Swift_CG2_LoadFlightData()


    % Inertial Data
    FD.I.g = 9.81;           % Gravity Constant
    FD.I.m = 598;           % Aircraft Mass (kg)
    FD.I.Ixx = 563;         % Aircraft Moments of Inertia (kg.m^2)
    FD.I.Iyy = 915;         % Aircraft Moments of Inertia (kg.m^2)
    FD.I.Izz = 1158;        % Aircraft Moments of Inertia (kg.m^2)
    FD.I.Ixz = 0;          % Aircraft Moments of Inertia (kg.m^2)
    
    % Geometric Data
    FD.Geo.S = 21.365;         % Platform Area (m^2)
    FD.Geo.c = 1.4142;         % Chord Length (m)
    FD.Geo.b = 8.4853;         % Wing Span (m)
    
    % Propeller Data
    FD.Prop.P_max = 97000; % Maximum engine power (Watts)
    FD.Prop.eta = 0.85; % Propeller efficiency
    
    % Control Data
    DtoR = pi/180;
    FD.CntrlLimit.Lower = [0;            % Throttle range (Fraction)
                                      -30*DtoR;     % Elevator range (rad)
                                      -30*DtoR;     % Aileron range (rad)
                                      -25*DtoR];    % Rudder range (rad)
    FD.CntrlLimit.Upper = [1;            % Throttle range (Fraction)
                                      30*DtoR;      % Elevator range (rad)
                                      30*DtoR;      % Aileron range (rad)
                                      25*DtoR];     % Rudder range (rad)

    % Aerodynamic Data (Reference CG: 22 % mac)
    FD.Aero.alpha_o = 0; 
    % Drag Coefficients
    FD.Aero.Cdo    =  0.039854;
    FD.Aero.k      =  0.060979;
    FD.Aero.CDdf   =  0;
    % Lift Coefficients
    FD.Aero.CLa  =  5.3053;
    FD.Aero.CLq  =  9.0957;
    FD.Aero.CLad = -1.4186;
    FD.Aero.CLde =  0.29198;
    FD.Aero.CLdf =  0;
    FD.Aero.CLo  = 0;
    % Side Force Coefficients
    FD.Aero.Cyb  = -0.26369;
    FD.Aero.Cybd = 0.1168; %%%%%
    FD.Aero.Cyp  = 0.060607; %%%%%
    FD.Aero.Cyr  =  0.069873;
    FD.Aero.Cyda =  -0.019444; %%%%%%
    FD.Aero.Cydr =  0.15934; %%%%%%
    % M Moment Coefficients
    FD.Aero.Cmo  =  0.03346;
    FD.Aero.Cma  = -0.341; %-0.641; original
    FD.Aero.Cmq  = -7.4702; % -9.4702; original
    FD.Aero.Cmad = 0.9646; %%%%%%%
    FD.Aero.Cmde = -0.717;
    FD.Aero.Cmdf = 0;
    % N Moment Coefficients
    FD.Aero.Cnb  =  0.018525; % 0.025525
    FD.Aero.Cnbd =  -0.0048;
    FD.Aero.Cnp  = -0.046552;
    FD.Aero.Cnr  = -0.058342; % -0.078342
    FD.Aero.Cnda =  -0.016169;
    FD.Aero.Cndr = -0.075229; %%%%%%
    % L Moment Coefficients
    FD.Aero.Clb  = -0.060463;
    FD.Aero.Clbd = 0.0068;
    FD.Aero.Clp  = -0.41978;
    FD.Aero.Clr  =  0.11668;
    FD.Aero.Clda = -0.32567; 
    FD.Aero.Cldr =  -0.0093965;

end      
