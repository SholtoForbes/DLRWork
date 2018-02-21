function [altdot,londot,latdot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,heating_rate] = SpaceLinerVehicleModel(states,controls,auxdata)

interp = auxdata.interp;

alt = states.alt;
lon = states.lon;
lat = states.lat;
v = states.v;
gamma = states.gamma;
zeta = states.zeta;

aoa = controls.Alpha;
bank = controls.eta;

% if isnan(alt)
%    alt = [1000; 60000]; 
% end

% =======================================================
% Vehicle Model
% =======================================================

A = auxdata.Stage2.A; % reference area (m^2)

%Gravity
g = 9.81;

% Change mass dependent on stage. 

m = auxdata.Stage2.mStruct+auxdata.Stage2.mFuel;


%======================================================

%% Flow =============================================================
%% Flow =============================================================
c = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
c(alt>=84000) = ppval(interp.c_spline,84000*ones(1,length(alt(alt>=84000))));
% c = c';

mach = v./c;
mach(alt>=84000) = v(alt>=84000)./c(alt>=84000);
% mach = mach';

rho = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data
% rho(alt>=84000) = ppval(interp.rho_spline,84000);
rho(alt>=84000) = ppval(interp.rho_spline,84000).*gaussmf (alt(alt>=84000), [10000 84000]);
% rho = rho';

% q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

% M = v./c; % Calculating Mach No (Descaled)
% 
% T0 = ppval(interp.T0_spline, alt); 
% 
% P0 = ppval(interp.P0_spline, alt);


% c = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
% mach = v./c;
% rho = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data

%% Aerodynamics
% interpolate coefficients


Cd = auxdata.interp.Stage2.Cd_spline(mach,rad2deg(aoa));
Cl = auxdata.interp.Stage2.Cl_spline(mach,rad2deg(aoa));   


%%%% Compute the drag and lift:

D = 0.5*Cd.*A.*rho.*v.^2;
L = 0.5*Cl.*A.*rho.*v.^2;

%% Thrust 

% rho_SL = ppval(interp.rho_spline,0);
% 
% 
% T1 = rho./rho_SL.*Stage1.T_SL + (1-rho./rho_SL).*Stage1.T_vac; % Thrust from stage 1
% T2 = 0;
% 
% T = T1 + T2;
% 
% Isp1 = rho./rho_SL.*Stage1.Isp_SL + (1-rho./rho_SL).*Stage1.Isp_vac;
% Isp2 = 0; 
% 
% Fueldt = T2./Isp2/g;

T=0;
Fueldt =0;


%%heating---------------------------
% From Conceptual Shape Optimization of Entry Vehicles, Dirkx & Mooj & NASA
% lecture

%using hot wall correction

% kappa = 1.7415e-4; % cited as sutton-graves, from nasa lecture
% Rn = 0.205; %effective nose radius (m) 
% 
% heating_rate = kappa*sqrt(rho./Rn).*v.^3; %W/m^2

%Heating model used in Tosca

R_N = 0.205; %effective nose radius (m) 

C = 20254.4;
rho_r = 1.225;
v_r = 10000;
R_Nr = 1;

heating_rate = C*sqrt(rho/rho_r*R_Nr/R_N).*(v/v_r).^3.05*1e4;


% 
% Q = zeros(1,length(time));
% Q(1) = 0;
% 
% for i = 1:length(dt_array)
%     Q(i+1) = heating_rate(i)*dt_array(i) + Q(i);
% end

%Motion in Geodetic Rotational Coordinates =================================================

[altdot,londot,latdot,gammadot,a,zetadot] = RotCoords(alt+auxdata.Re,lon,lat,gamma,v,zeta,L,D,T,m,aoa,bank);

% Aero Data =============================================================

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

v_H = v.*cos(gamma);

% =========================================================================
end








