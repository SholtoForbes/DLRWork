function [altdot,xidot,phidot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,T,Isp1,Isp2,m,heating_rate,total_acceleration,Cl,Cd] = SpaceLinerVehicleModel(t,phase,throttle,auxdata,stage)

% STF = 0.6; %Staging Time Fraction


alt     = phase.state(:,1);
lon     = phase.state(:,2);
lat     = phase.state(:,3);
v       = phase.state(:,4);
gamma   = phase.state(:,5);
zeta    = phase.state(:,6);
mFuel   = phase.state(:,7);
Alpha = phase.state(:,8);
eta = phase.state(:,9);

interp = auxdata.interp;




Stage1 = auxdata.Stage1;
Stage2 = auxdata.Stage2;


% if isnan(timeF)
%    t =  0:1/(length(alt)-1):1;
%    timeF = 1;
% end

% if isnan(alt)
%    alt = [1000; 60000]; 
% end


% Stage = 1;

% =======================================================
% Vehicle Model
% =======================================================

A1 = auxdata.Stage1.A; % reference area (m^2)
A2 = auxdata.Stage2.A; % reference area (m^2)
%Gravity
g = 9.81;

% Change mass dependent on stage. 
if stage == 1
    m = auxdata.Stage1.mStruct+mFuel+auxdata.Stage2.mStruct; 
elseif stage == 2
    m = auxdata.Stage2.mStruct+mFuel;
end

% m  = auxdata.Stage1.mStruct+mFuel +auxdata.Stage2.mStruct;
% 
% m  = auxdata.Stage2.mStruct+mFuel ;

%======================================================

%% Flow =============================================================
c = ppval(interp.c_spline,alt); % Calculate speed of sound using atmospheric data
c(alt>=84000) = ppval(interp.c_spline,84000*ones(1,length(alt(alt>=84000))));
% c = c';

mach = v./c;
mach(alt>=84000) = v(alt>=84000)./c(alt>=84000);
% mach = mach';

rho = ppval(interp.rho_spline,alt); % Calculate density using atmospheric data
% rho(alt>=84000) = ppval(interp.rho_spline,84000);
rho(alt>=84000) = ppval(interp.rho_spline,84000).*gaussmf (alt(alt>=84000), [1000 84000]);
% rho = rho';

% q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure

% M = v./c; % Calculating Mach No (Descaled)

% T0 = ppval(interp.T0_spline, alt); 
% P0 = ppval(interp.P0_spline, alt);




%% Aerodynamics
% interpolate coefficients

rho_SL = ppval(interp.rho_spline,0);

if stage == 1
Cd  = auxdata.interp.Stage1.Cd_spline(mach ,rad2deg(Alpha ));
Cl  = auxdata.interp.Stage1.Cl_spline(mach ,rad2deg(Alpha ));
D  = 0.5*Cd.*A1.*rho.*v.^2;
L  = 0.5*Cl.*A1.*rho.*v.^2;
  T1  = rho ./rho_SL.*Stage1.T_SL + (1-rho ./rho_SL).*Stage1.T_vac; % Thrust from stage 1
    T2  = rho ./rho_SL.*Stage2.T_SL + (1-rho ./rho_SL).*Stage2.T_vac;
    

    
    Isp1  = rho ./rho_SL.*Stage1.Isp_SL + (1-rho ./rho_SL).*Stage1.Isp_vac;
    Isp2  = rho ./rho_SL.*Stage2.Isp_SL + (1-rho ./rho_SL).*Stage2.Isp_vac;
  Fueldt  = (T1 ./Isp1 /g + T2 ./Isp2 /g).*throttle;


elseif stage == 2
    Cd  = auxdata.interp.Stage2.Cd_spline(mach ,rad2deg(Alpha ));
Cl  = auxdata.interp.Stage2.Cl_spline(mach ,rad2deg(Alpha ));
D  = 0.5*Cd.*A2.*rho.*v.^2;
L  = 0.5*Cl.*A2.*rho.*v.^2;

  
    T1  = 0; % Thrust from stage 1
    T2  = rho ./rho_SL.*Stage2.T_SL + (1-rho ./rho_SL).*Stage2.T_vac;
    

    Isp1  = 0;
    Isp2  = rho ./rho_SL.*Stage2.Isp_SL + (1-rho ./rho_SL).*Stage2.Isp_vac; 
     Fueldt  = (T2 ./Isp2 /g).*throttle;


end




%% Thrust 

     
    T1 = T1.*throttle;
    T2 = T2.*throttle;
    
    T = T1 + T2;

%Motion in Geodetic Rotational Coordinates =================================================
% eta = 0;
[altdot,xidot,phidot,gammadot,a,zetadot,total_lift] = RotCoords(alt'+auxdata.Re,lon',lat',gamma',v',zeta',L',D',T',m',Alpha',eta');

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

% total_acceleration = a;

total_acceleration = sqrt(a.^2 + (total_lift./m').^2); % Lateral and axial acceleration
end








