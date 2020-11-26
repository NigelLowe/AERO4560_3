%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Non linear integration
% Function to perform simple Euler integration of state equations 
% for generic flight dynamic model
%
% (c) Peter W. Gibbens & S. Dumble, 1 March, 2011.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [X_out] = aero4560_eulerNL(DT, X, Xg, U, FlightData)
    % change folder to previous to get .p functions
    Current_Folder = pwd;
    TopFolder = fileparts(pwd);
    cd(TopFolder);


    % Obtain state derivatives (Assume Xdot = 0)
    [ForceCoeff, MomentCoeff] = aero4560_aero(X, Xg, zeros(12,1), U, FlightData);
    [Xdot]  = aero4560_motion(X, ForceCoeff, MomentCoeff, FlightData);
    % Obtain state derivatives (Iterate Once so that Xdot Converges)
    [ForceCoeff, MomentCoeff] = aero4560_aero(X, Xg, Xdot, U, FlightData);
    [Xdot]  = aero4560_motion(X, ForceCoeff, MomentCoeff, FlightData);

    % Integrate states
    X_out = X + Xdot*DT;
    
    % move back to original folder
    cd(Current_Folder);
end


