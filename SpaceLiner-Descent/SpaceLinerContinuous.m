function phaseout = SpaceLinerContinuous(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics
states.alt     = input.phase(1).state(:,1);
states.lon     = input.phase(1).state(:,2);
states.lat     = input.phase(1).state(:,3);
states.v       = input.phase(1).state(:,4);
states.gamma   = input.phase(1).state(:,5);
states.zeta    = input.phase(1).state(:,6);

controls.Alpha = input.phase(1).state(:,7); % Note the 'controls' class here is not the same as the 'control' defined for GPOPS (these are the vehicle controls, not the control theory controls)
% controls.eta   = input.phase(1).state(:,8);

controls.throttle  = 1;

Alphadot1  = input.phase(1).control(:,1);
etadot1 = input.phase(1).control(:,2);

time = input.phase(1).time;

auxdata = input.auxdata;

[altdot,londot,latdot,gammadot,vdot,zetadot, q1, M, D, rho,L,Fueldt,heating_rate] = SpaceLinerVehicleModel(states,controls,auxdata);

phaseout(1).dynamics  = [altdot, londot, latdot, vdot, gammadot, zetadot, Alphadot1];


total_acceleration = sqrt(vdot.^2 + (states.v.*gammadot).^2 + (states.v.*zetadot).^2)/9.81;

    
% phaseout(1).path = [q1];

phaseout(1).path = [q1, heating_rate, total_acceleration];



states.lon(isnan(states.lon)) = 0;
states.lat(isnan(states.lat)) = 0;

states.lon = states.lon + auxdata.lon0;
states.lon(states.lon > pi) = states.lon(states.lon > pi) - 2*pi;
states.lon(states.lon < -pi) = states.lon(states.lon < -pi) + 2*pi;

% isLand = auxdata.LandSpline(rad2deg(states.lat),rad2deg(states.lon));
% 
% isLand(rad2deg(states.lat) < -60) = 0; % Set antarctica to not be relevant
% 
% AltOverLandCost = zeros(length(states.alt),1);
% AltOverLandCost(isLand==1) = 80000-states.alt(isLand==1); % Find altitudes of every point over land, maximise
% AltOverLandCost(states.alt>80000) = 0;
% % The 'target' must be greater than the max alt, or go to 0 above target,
% % to prevent it being better to fly high over land, than over water
% 
% phaseout(1).integrand = AltOverLandCost.^2/1e6;
% % phaseout(1).integrand = isLand;


AltCost = (80000-states.alt)/80000;

pop = auxdata.PopInterp(rad2deg(states.lon),rad2deg(states.lat));

% popCost = pop.*AltCost; % for flights which go over large amounts of
% population

popCost = pop;

% phaseout(1).integrand = pop;

% Scaling of heating must change if different population density acfuracy
% is used, because population density is per cell, which will change in magnitude.
% phaseout(1).integrand = pop + heating_rate/1e6;
% phaseout(1).integrand = pop + heating_rate/1e5;
phaseout(1).integrand = popCost + heating_rate/1e5;
end

%======================================================