function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,Isp1,Isp2,q1,flap_deflection] = SpaceLinerVehicleModel(states,controls,auxdata,Stage)

interp = auxdata.interp;

alt = states.alt;
lon = states.lon;
lat = states.lat;
v = states.v;
gamma = states.gamma;
zeta = states.zeta;
mFuel = states.mFuel;

Alpha = controls.Alpha;
eta = controls.eta;


% =======================================================
% Vehicle Model
% =======================================================

A = auxdata.A; % reference area (m^2)

%Gravity
g = 9.81;

% Change mass dependent on stage. 
if Stage == 1
    m = auxdata.Stage1.mStruct+mFuel+auxdata.Stage2.mStruct; 
else
    m = auxdata.Stage2.mStruct+mFuel;
end

%======================================================

%% Flow =============================================================
c = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
mach = v./c;
rho = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

T0 = ppval(interp.T0_spline, alt); 

P0 = ppval(interp.P0_spline, alt);

%% Aerodynamics
% interpolate coefficients

if Stage == 1
    Cd = auxdata.interp.Stage1.Cd_spline(mach,rad2deg(alpha));
    Cl = auxdata.interp.Stage1.Cl_spline(mach,rad2deg(alpha));
else
    Cd = auxdata.interp.Stage2.Cd_spline(mach,rad2deg(alpha));
    Cl = auxdata.interp.Stage2.Cl_spline(mach,rad2deg(alpha));   
end

%%%% Compute the drag and lift:

D = 0.5*Cd.*A.*rho.*v.^2;
L = 0.5*Cl.*A.*rho.*v.^2;

%% Thrust 

rho_SL = ppval(interp.rho_spline,0);

if Stage == 1
    T1 = rho./rho_SL.*Stage1.T_SL + (1-rho./rho_SL).*Stage1.T_vac; % Thrust from stage 1
    T2 = rho./rho_SL.*Stage2.T_SL + (1-rho./rho_SL).*Stage2.T_vac;
    
    T = T1 + T2;
    
    Isp1 = rho./rho_SL.*Stage1.Isp_SL + (1-rho./rho_SL).*Stage1.Isp_vac;
    Isp2 = rho./rho_SL.*Stage2.Isp_SL + (1-rho./rho_SL).*Stage2.Isp_vac;

    Fueldt = T1./Isp1/g + T2./Isp2/g;
else
    T1 = rho./rho_SL.*Stage1.T_SL + (1-rho./rho_SL).*Stage1.T_vac; % Thrust from stage 1
    T2 = 0;
    
    T = T1 + T2;
    
    Isp1 = rho./rho_SL.*Stage1.Isp_SL + (1-rho./rho_SL).*Stage1.Isp_vac;
    Isp2 = 0; 

    Fueldt = T2./Isp2/g;
end

%Motion in Geodetic Rotational Coordinates =================================================

[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoordsReturn(alt+auxdata.Re,xi,phi,gamma,v,zeta,L,D,T,m,alpha,eta);

% Aero Data =============================================================

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

v_H = v.*cos(gamma);

% =========================================================================
end








