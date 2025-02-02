function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,Isp1,Isp2,m,heating_rate] = SpaceLinerVehicleModel(t,states,controls,throttle,auxdata,timeF)

STF = 0.6; %Staging Time Fraction



interp = auxdata.interp;

alt = states.alt;
lon = states.lon;
lat = states.lat;
v = states.v;
gamma = states.gamma;
zeta = states.zeta;
mFuel = states.mFuel;

Alpha = controls.Alpha;
% eta = controls.eta;

Stage1 = auxdata.Stage1;
Stage2 = auxdata.Stage2;


if isnan(timeF)
   t =  0:1/(length(alt)-1):1;
   timeF = 1;
end

if isnan(alt)
   alt = [1000; 60000]; 
end


Stage = 1;

% =======================================================
% Vehicle Model
% =======================================================

A1 = auxdata.Stage1.A; % reference area (m^2)
A2 = auxdata.Stage2.A; % reference area (m^2)
%Gravity
g = 9.81;

% Change mass dependent on stage. 
% if Stage == 1
%     m = auxdata.Stage1.mStruct+mFuel+auxdata.Stage2.mStruct; 
% else
%     m = auxdata.Stage2.mStruct+mFuel;
% end

m(t<timeF*STF) = auxdata.Stage1.mStruct+mFuel(t<timeF*STF)+auxdata.Stage2.mStruct; 
m(t>=timeF*STF) = auxdata.Stage2.mStruct+mFuel(t>=timeF*STF);

%======================================================

%% Flow =============================================================
c(alt<84000) = ppval(interp.c_spline,alt(alt<84000)); % Calculate speed of sound using atmospheric data
c(alt>=84000) = ppval(interp.c_spline,84000*ones(1,length(alt(alt>=84000))));
c = c';

mach(alt<84000) = v(alt<84000)./c(alt<84000);
mach(alt>=84000) = v(alt>=84000)./c(alt>=84000);
mach = mach';

rho(alt<84000) = ppval(interp.rho_spline,alt(alt<84000)); % Calculate density using atmospheric data
rho(alt>=84000) = ppval(interp.rho_spline,84000);
rho = rho';

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

% M = v./c; % Calculating Mach No (Descaled)

% T0 = ppval(interp.T0_spline, alt); 
% P0 = ppval(interp.P0_spline, alt);




%% Aerodynamics
% interpolate coefficients

% if Stage == 1
%     Cd = auxdata.interp.Stage1.Cd_spline(mach,rad2deg(Alpha));
%     Cl = auxdata.interp.Stage1.Cl_spline(mach,rad2deg(Alpha));
% else
%     Cd = auxdata.interp.Stage2.Cd_spline(mach,rad2deg(Alpha));
%     Cl = auxdata.interp.Stage2.Cl_spline(mach,rad2deg(Alpha));   


Cd(t<timeF*STF) = auxdata.interp.Stage1.Cd_spline(mach(t<timeF*STF),rad2deg(Alpha(t<timeF*STF)));
Cd(t>=timeF*STF) = auxdata.interp.Stage2.Cd_spline(mach(t>=timeF*STF),rad2deg(Alpha(t>=timeF*STF)));

Cl(t<timeF*STF) = auxdata.interp.Stage1.Cl_spline(mach(t<timeF*STF),rad2deg(Alpha(t<timeF*STF)));
Cl(t>=timeF*STF) = auxdata.interp.Stage2.Cl_spline(mach(t>=timeF*STF),rad2deg(Alpha(t>=timeF*STF)));

%%%% Compute the drag and lift:
D(t<timeF*STF) = 0.5*Cd(t<timeF*STF)'.*A1.*rho(t<timeF*STF).*v(t<timeF*STF).^2;
D(t>=timeF*STF) = 0.5*Cd(t>=timeF*STF)'.*A2.*rho(t>=timeF*STF).*v(t>=timeF*STF).^2;
D = D';

L(t<timeF*STF) = 0.5*Cl(t<timeF*STF)'.*A1.*rho(t<timeF*STF).*v(t<timeF*STF).^2;
L(t>=timeF*STF) = 0.5*Cl(t>=timeF*STF)'.*A2.*rho(t>=timeF*STF).*v(t>=timeF*STF).^2;
L = L';

%% Thrust 

rho_SL = ppval(interp.rho_spline,0);

% if Stage == 1
%     T1 = rho./rho_SL.*Stage1.T_SL + (1-rho./rho_SL).*Stage1.T_vac; % Thrust from stage 1
%     T2 = rho./rho_SL.*Stage2.T_SL + (1-rho./rho_SL).*Stage2.T_vac;
%     
%     T = (T1 + T2).*throttle;
%     
%     Isp1 = rho./rho_SL.*Stage1.Isp_SL + (1-rho./rho_SL).*Stage1.Isp_vac;
%     Isp2 = rho./rho_SL.*Stage2.Isp_SL + (1-rho./rho_SL).*Stage2.Isp_vac;
% 
%     Fueldt = (T1./Isp1/g + T2./Isp2/g).*throttle;
%     
%     
% else
%     T1 = rho./rho_SL.*Stage1.T_SL + (1-rho./rho_SL).*Stage1.T_vac; % Thrust from stage 1
%     T2 = 0;
%     
%     T = T1 + T2;
%     
%     Isp1 = rho./rho_SL.*Stage1.Isp_SL + (1-rho./rho_SL).*Stage1.Isp_vac;
%     Isp2 = 0; 
% 
%     Fueldt = T2./Isp2/g;
% end


    T1(t<timeF*STF) = rho(t<timeF*STF)./rho_SL.*Stage1.T_SL + (1-rho(t<timeF*STF)./rho_SL).*Stage1.T_vac; % Thrust from stage 1
    T2(t<timeF*STF) = rho(t<timeF*STF)./rho_SL.*Stage2.T_SL + (1-rho(t<timeF*STF)./rho_SL).*Stage2.T_vac;
    

    
    Isp1(t<timeF*STF) = rho(t<timeF*STF)./rho_SL.*Stage1.Isp_SL + (1-rho(t<timeF*STF)./rho_SL).*Stage1.Isp_vac;
    Isp2(t<timeF*STF) = rho(t<timeF*STF)./rho_SL.*Stage2.Isp_SL + (1-rho(t<timeF*STF)./rho_SL).*Stage2.Isp_vac;


    T1(t>=timeF*STF) = 0; % Thrust from stage 1
    T2(t>=timeF*STF) = rho(t>=timeF*STF)./rho_SL.*Stage2.T_SL + (1-rho(t>=timeF*STF)./rho_SL).*Stage2.T_vac;
    

    Isp1(t>=timeF*STF) = 0;
    Isp2(t>=timeF*STF) = rho(t>=timeF*STF)./rho_SL.*Stage2.Isp_SL + (1-rho(t>=timeF*STF)./rho_SL).*Stage2.Isp_vac; 
    

    Fueldt(t<timeF*STF) = (T1(t<timeF*STF)./Isp1(t<timeF*STF)/g + T2(t<timeF*STF)./Isp2(t<timeF*STF)/g).*throttle(t<timeF*STF)';
    Fueldt(t>=timeF*STF) = (T2(t>=timeF*STF)./Isp2(t>=timeF*STF)/g).*throttle(t>=timeF*STF)';

    T1 = T1.*throttle';
    T2 = T2.*throttle';
    
    T = T1 + T2;

%Motion in Geodetic Rotational Coordinates =================================================
eta = 0;
[altdot,xidot,phidot,gammadot,a,zetadot] = RotCoords(alt+auxdata.Re,lon,lat,gamma,v,zeta,L,D,T',m',Alpha,eta);

% Aero Data =============================================================

q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

M = v./c; % Calculating Mach No (Descaled)

v_H = v.*cos(gamma);

%Heating model used in Tosca

R_N = 0.205; %effective nose radius (m) 

C = 20254.4;
rho_r = 1.225;
v_r = 10000;
R_Nr = 1;

heating_rate = C*sqrt(rho/rho_r*R_Nr/R_N).*(v/v_r).^3.05*1e4;
% =========================================================================
end








