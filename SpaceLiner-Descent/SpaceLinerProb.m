%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rocket Flight Optimiser
% By Sholto Forbes-Spyratos
% Utilises the GPOPS-II proprietary optimisation software
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Timestamp = datestr(now,30)

%% Atmosphere Data %%======================================================
% Fetch atmospheric data and compute interpolation splines.

Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
auxdata.interp.rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data
auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data
auxdata.interp.T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 
auxdata.interp.P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 

%% Import Vehicle Config Data %%============================

Stage2.A = 461; %Reference Area in m² 


Stage2.mStruct = 134361.1;
% Stage2.mFuel = 5384.7;

% Japan-germany
Stage2.mFuel = 172586-Stage2.mStruct;

Stage2.T_SL = 1830*2*1e3; %for mixture ratio 6 with all boosters active
Stage2.Isp_SL = 363;

Stage2.T_vac = 2268*2*1e3;
Stage2.Isp_vac = 449;

auxdata.Stage2 = Stage2;


%%
auxdata.Re   = 6371203.92;                     % Equatorial Radius of Earth (m)

%%
% =========================================================================
% SET RUN MODE
% =========================================================================
% Change const to set the target of the simulation. Much of the problem
% definition changes with const.

% const = 1x: No end constraint, used for optimal trajectory calculation
% const = 1: 50kPa limit, 12: 55 kPa limit, 13: 45 kPa limit, 14: 50kPa limit & 10% additional drag

const = 1
auxdata.const = const;

%% Aerodynamic Data
% Fetch aerodynamic data and compute interpolation splines.

Stage2_MList = [0.7  0.9  1.1  2.0  4.0  6.0  10.0  14.0  17.9  18.1  22.0  26.0];
Stage2_aoaList = [0.00   1.00   2.00   3.00   4.00   5.00   6.00   7.00   8.00   9.00  10.00  11.00  12.00  13.00  14.00  15.00  16.00  17.00  18.00  19.00  20.00  21.00  22.00  23.00  24.00  25.00  26.00  27.00  28.00  29.00  30.00];
Stage2_Cl = importdata('aero_Orbiter_Lift');
Stage2_Cd = importdata('aero_Orbiter_Drag');

Stage2_MList_Full = zeros(length(Stage2_Cl),1);
Stage2_aoaList_Full = zeros(length(Stage2_Cl),1);

for i = 1:length(Stage2_Cl)
    Stage2_aoaList_Full(i) = Stage2_aoaList(i - length(Stage2_aoaList)*floor((i-1)/length(Stage2_aoaList))); 
    Stage2_MList_Full(i) = Stage2_MList(floor((i-1)/length(Stage2_aoaList)) + 1);
end

[Stage2_MGrid,Stage2_aoaGrid] = meshgrid(Stage2_MList,Stage2_aoaList);

Stage2_Cl_scattered = scatteredInterpolant(Stage2_MList_Full,Stage2_aoaList_Full,Stage2_Cl);
Stage2_Cd_scattered = scatteredInterpolant(Stage2_MList_Full,Stage2_aoaList_Full,Stage2_Cd);

for i = 1:length(Stage2_MList)
    for j = 1:length(Stage2_aoaList)
        Stage2_Cl_Grid(j,i) = Stage2_Cl_scattered(Stage2_MGrid(j,i),Stage2_aoaGrid(j,i));
        Stage2_Cd_Grid(j,i) = Stage2_Cd_scattered(Stage2_MGrid(j,i),Stage2_aoaGrid(j,i));
    end
end

%% Create Gridded Interpolants

auxdata.interp.Stage2.Cl_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cl_Grid','spline');
auxdata.interp.Stage2.Cd_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cd_Grid','spline');


%% Get land interpolants to determine if over land
load LandSpline

auxdata.LandSpline = LandSpline;


%% Get population interpolation
load PopInterp

auxdata.PopInterp = PopInterp;

%% Set Bounds %%========================================================

%% If initial heading angle is bounded

zeta0 = deg2rad(70.18);

%%

altMin = 1;
altMax = 90000;


vMin = 1;
vMax = 10000;

gammaMin = -deg2rad(80);
gammaMax = deg2rad(80);

zetaMin = -2*pi;
% zetaMin = 0;
zetaMax = 2*pi;

lonMin = -2*pi;         
lonMax = 2*pi;

% latMin = -pi/2+0.0000001;  
latMin = 0;  %for japan-germany actual

latMax = pi/2-0.0000001;
% latMax = pi;

aoaMin = 0;  
aoaMax = deg2rad(25);

bankMin = deg2rad(-50); 
bankMax =   deg2rad(50);

% bankMin = deg2rad(-60); 
% bankMax =   deg2rad(60);

% Initial Conditions

% lat0 = deg2rad(-23.3791); % Rockhampton
% lon0 = deg2rad(150.5100); % Rockhampton

% lat0 = deg2rad(30.373796); % South Japan
% lon0 = deg2rad(130.95852);

% lat0 = deg2rad(37.533151); % Close to Suzu, middle north Japan
% lon0 = deg2rad(137.276039);

% lat0 = deg2rad(38.849086); % Taiyuan Satellite Launch Center
% lon0 = deg2rad(111.608497);

% lat0 = deg2rad(40.9675); % Jiuquan Satellite Launch Center (Szechuan)
% lon0 = deg2rad(100.278611); 

% lat0 = deg2rad(53.77); % Germany
% lon0 = deg2rad(8.6359);% Germany

% lat0 = deg2rad(28.469159); %Cape Canaveral
% lon0 = deg2rad(-80.509981);

% lat0 = deg2rad(45.522630); %Cape Soya, Hokkaido, northernmost tip of Japan
% lon0 = deg2rad(141.937105);

% lat0 = deg2rad(38.745095); %north korea
% lon0 = deg2rad(128.272375);

% lat0 = deg2rad(37.235705); %south korea
% lon0 = deg2rad(129.356220);

% lat0 = deg2rad(-32.711986); % Close to cape town 
% lon0 = deg2rad(17.927735);

% lat0 = deg2rad(35.49196); % South of south korea
% lon0 = deg2rad(129.444331);


% lat0 = deg2rad(35.528217); %South Japan
% lon0 = deg2rad(133.569646);

lat0 = deg2rad(45.37); %South Japan after true launch
lon0 = deg2rad(137.71);

auxdata.lon0 = lon0;

% alt0 = 70000;
% v0 = 7000;
% gamma0 = 0;

%japan-germany
alt0 = 73235;
v0 = 6523;
gamma0 = deg2rad(0.128);

% End conditions
altFMin = 100;
altFMax = 1000;

latF = deg2rad(53.77); % Germany
lonF = deg2rad(8.6359);% Germany

% latF = deg2rad(45.522630); %Cape Soya, Hokkaido, northernmost tip of Japan
% lonF = deg2rad(141.937105);


% latF = deg2rad(28.469159); %Cape Canaveral
% lonF = deg2rad(-80.509981);

% latF = deg2rad(-23.962496); % Brazil
% lonF = deg2rad(-45.237146);

% latF = deg2rad(40); %Test
% lonF = deg2rad(-30);

% latF = deg2rad(-23.3791); % Rockhampton
%  lonF = deg2rad(150.5100); % Rockhampton

% latF = deg2rad(37.533151); % Close to Suzu, middle north Japan
% lonF = deg2rad(137.276039);

% latF = deg2rad(-32.711986); % Close to cape town 
% lonF = deg2rad(17.927735);

% lonF = deg2rad(218);
% Primal Bounds
bounds.phase(1).state.lower = [altMin, lonMin, latMin, vMin, gammaMin, zetaMin, aoaMin, bankMin];
bounds.phase(1).state.upper = [altMax, lonMax, latMax, vMax, gammaMax, zetaMax, aoaMax, bankMax];

% Initial States
% bounds.phase(1).initialstate.lower = [alt0,0, lat0, v0, gamma0, zetaMin, aoaMin, bankMin] ;
% bounds.phase(1).initialstate.upper = [alt0,0, lat0, v0, gamma0, zetaMax, aoaMax, bankMax];

% bounds.phase(1).initialstate.lower = [alt0,0, lat0, v0-1000, gamma0, zetaMin, aoaMin, bankMin] ;
% bounds.phase(1).initialstate.upper = [alt0,0, lat0, v0, gamma0, zetaMax, aoaMax, bankMax];

bounds.phase(1).initialstate.lower = [alt0,0, lat0, v0, gamma0, zeta0, aoaMin, bankMin] ;
bounds.phase(1).initialstate.upper = [alt0,0, lat0, v0, gamma0, zeta0, aoaMax, bankMax];


% End States
%Starting East heading east
% For aus-germany, japan-germany
bounds.phase(1).finalstate.lower = [altFMin, lonF-lon0+2*pi, latF, vMin, deg2rad(-20), zetaMin, aoaMin, 0];
bounds.phase(1).finalstate.upper = [altFMax, lonF-lon0+2*pi, latF, 150, deg2rad(20), zetaMax, aoaMax, 0];

%Starting West heading west
% for germany-japan
% bounds.phase(1).finalstate.lower = [altF, lonF-lon0-2*pi, latF, vMin, -0.05, zetaMin, aoaMin, 0];
% bounds.phase(1).finalstate.upper = [altF, lonF-lon0-2*pi, latF, 50, 0.05, zetaMax, aoaMax, 0];

% If not going outside of degree bounds
% for japan-florida, florida-japan
% bounds.phase(1).finalstate.lower = [altF, lonF-lon0, latF, vMin, -0.05, zetaMin, aoaMin, 0];
% bounds.phase(1).finalstate.upper = [altF, lonF-lon0, latF, 50, 0.05, zetaMax, aoaMax, 0];

% Control Bounds
bounds.phase(1).control.lower = [deg2rad(-.1), deg2rad(-1)];
bounds.phase(1).control.upper = [deg2rad(.1), deg2rad(1)];

% Time Bounds
bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower =0;
bounds.phase(1).finaltime.upper = 10000;

%Integral Bounds
% bounds.phase.integral.lower = 0;
% bounds.phase.integral.upper = 100000;

bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 1e9;

%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.
% 
% bounds.phase(1).path.lower = [0, -1];
% bounds.phase(1).path.upper = [60000, 100];

% bounds.phase(1).path.lower = [0, 0, -2.5];
% bounds.phase(1).path.upper = [60000, 2e6, 2.5];

bounds.phase(1).path.lower = [0, 0, -2.5]; 
bounds.phase(1).path.upper = [60000, 1.5e6, 2.5];

%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(1).state(:,1)   = [alt0;alt0];
guess.phase(1).state(:,2)   = [0;lonF-lon0+2*pi]; 
% guess.phase(1).state(:,2)   = [0;lonF-lon0-2*pi];
% guess.phase(1).state(:,2)   = [0;lonF-lon0];
guess.phase(1).state(:,3)   = [lat0;latF];
% guess.phase(1).state(:,3)   = [lat0;latF-0.5];
guess.phase(1).state(:,4)   = [v0;vMin];

% This is gamma
guess.phase(1).state(:,5)   = [0;0]; %

% guess.phase(1).state(:,6)   = [deg2rad(180);deg2rad(180)];
% guess.phase(1).state(:,6)   = [0;0]; % For eastward
% guess.phase(1).state(:,6)   = [deg2rad(90);deg2rad(-45)]; % Aus-Germany
% guess.phase(1).state(:,6)   = [deg2rad(60);deg2rad(-40)]; % Japan-Germany
% guess.phase(1).state(:,6)   = [deg2rad(-10);deg2rad(90)]; %Aus-Canaveral
% guess.phase(1).state(:,6)   = [deg2rad(110);deg2rad(270)]; %Japan-Florida
% guess.phase(1).state(:,6)   = [deg2rad(90);deg2rad(270)]; %Aus-Canaveral west test
% guess.phase(1).state(:,6)   = [deg2rad(70);deg2rad(-80)];%Canaveral-Japan
% guess.phase(1).state(:,6)   = [deg2rad(100);deg2rad(180)]; % Germany-Japan
% guess.phase(1).state(:,6)   = [deg2rad(-10);deg2rad(-10)]; %Canaveral - cape town
% guess.phase(1).state(:,6)   = [deg2rad(110);deg2rad(110)]; %cape town-canaveral
% guess.phase(1).state(:,6)   = [deg2rad(-90);deg2rad(100)]; % Aus-Brazil
% guess.phase(1).state(:,6)   = [deg2rad(80);deg2rad(-60)]; % Korea-Germany
guess.phase(1).state(:,6)   = [deg2rad(70);deg2rad(-60)]; % japan-Germany actual


 guess.phase(1).state(:,7)   = [10*pi/180; 10*pi/180];
guess.phase(1).state(:,8)   = [0;0];
% guess.phase(1).state(:,8)   = [deg2rad(50);deg2rad(50)];

guess.phase(1).control      = [[0;0],[0;0]];
guess.phase(1).time          = [0;4000];

guess.phase.integral = 0;

% Tire stages together
% bounds.eventgroup(1).lower = [zeros(1,10)];
% bounds.eventgroup(1).upper = [zeros(1,10)]; 


%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
% mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 5;
mesh.colpointsmin = 3;
mesh.colpointsmax = 200;
mesh.tolerance    = 1e-5;


%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'Reusable-Launch-Vehicle-Entry-Problem';
setup.functions.continuous           = @SpaceLinerContinuous;
setup.functions.endpoint             = @SpaceLinerEndpoint;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.maxiterations = 1000;
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
% setup.scales.method                  = 'automatic-guessUpdate';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%


output = gpops2(setup);

%%

EndTime = datestr(now,30) % Display the ending time

% =========================================================================
% Assign the primal variables
alt = output.result.solution.phase(1).state(:,1).';
lon = output.result.solution.phase(1).state(:,2).';
lat = output.result.solution.phase(1).state(:,3).';
v = output.result.solution.phase(1).state(:,4).'; 
gamma = output.result.solution.phase(1).state(:,5).'; 
zeta = output.result.solution.phase(1).state(:,6).';
Alpha = output.result.solution.phase(1).state(:,7).';
eta = output.result.solution.phase(1).state(:,8).';


omegadot  = output.result.solution.phase(1).control.'; 


time = output.result.solution.phase(1).time.';



states.alt     = output.result.solution.phase(1).state(:,1);
states.lon     = output.result.solution.phase(1).state(:,2);
states.lat     = output.result.solution.phase(1).state(:,3);
states.v       = output.result.solution.phase(1).state(:,4);
states.gamma   = output.result.solution.phase(1).state(:,5);
states.zeta    = output.result.solution.phase(1).state(:,6);

controls.Alpha = output.result.solution.phase(1).state(:,7); % Note the 'controls' class here is not the same as the 'control' defined for GPOPS (these are the vehicle controls, not the control theory controls)
controls.eta   = output.result.solution.phase(1).state(:,8);

[altdot,londot,latdot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,heating_rate] = SpaceLinerVehicleModel(states,controls,auxdata);

total_acceleration = sqrt(a.^2 + (states.v.*gammadot).^2 + (states.v.*zetadot).^2)/9.81;

pop = auxdata.PopInterp(rad2deg(states.lon),rad2deg(states.lat));

figure(201)
subplot(10,1,1)
hold on
plot(time,alt)
title('alt')
subplot(10,1,2)
hold on
plot(time,v)
title('v')

subplot(10,1,3)
hold on
plot(time,lon)
title('lon')

subplot(10,1,4)

hold on
plot(time,lat)
title('lat')

subplot(10,1,5)
hold on
plot(time,rad2deg(gamma))
title('gamma')

subplot(10,1,6)
hold on
plot(time,a/9.81)
title('acceleration (g)')

subplot(10,1,7)
hold on
plot(time,rad2deg(Alpha))
title('alpha')

subplot(10,1,8)
hold on
plot(time,rad2deg(eta))
title('eta')

subplot(10,1,9)
hold on
plot(time,heating_rate)
title('heat rate')

subplot(10,1,10)
hold on
plot(time,total_acceleration)
title('Acceleration, Experienced (g)')

figure(230)
hold on

axesm('pcarree','Origin',[0 rad2deg(lon0)+60 0])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
plotm(rad2deg(lat),rad2deg(lon)+rad2deg(lon0))

% plot3(lon-deg2rad(60),lat,alt/max(alt)) ; % Normalised altitude plot

cities = shaperead('worldcities', 'UseGeoCoords', true);
lats = extractfield(cities,'Lat');
lons = extractfield(cities,'Lon');
geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r',...
        'MarkerSize', 2)

figure(231)
% 
load topo
axesm ortho

axesm('globe','Geoid',earthRadius)
meshm(topo,topolegend);
demcmap(topo)

plot3m(rad2deg(states.lat),rad2deg(states.lon+lon0),'r','LineWidth',2)

geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'y',...
        'MarkerFaceColor', 'y',...
        'MarkerSize', 2)
% =========================================================================
% 
% 
% nodes = length(alt)
% 
% 
% 
% [~,~,~,~,~,~, q1, M1, Fd1, rho1,L1,Fueldt1,T1,Isp1,q11,flapdeflection] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,lat,lon,Alpha,eta,1, mFuel,mFuel(1),mFuel(end), 1);
% [~,~,~,~,~,~, q2, M2, Fd2, rho2,L2,Fueldt2,T2,Isp2,q12,flapdeflection2] = VehicleModelCombined(gamma2, alt2, v2,auxdata,zeta2,lat2,lon2,Alpha2,eta2,throttle2, mFuel2,0,0, 0);
% 
% throttle2(M2<5.0) = 0; % remove nonsense throttle points
% 
% % figure out horizontal motion
% H(1) = 0;
% for i = 1:nodes-1
% H(i+1) = v(i)*(time(i+1) - time(i))*cos(gamma(i)) + H(i);
% end
% 
% % Separation_LD = lift(end)/Fd(end)
% 
% figure(2010)
% 
% subplot(5,5,[1,10])
% hold on
% plot(H, alt)
% % plot(H(algorithm.nodes(1)), V(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
% title('Trajectory (m)')
% 
% dim = [.7 .52 .2 .2];
% annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Second Stage Fuel Used: ' num2str(1562 - mFuel(end)) ' kg']},'FitBoxToText','on');  
% 
% 
% subplot(5,5,11)
% hold on
% plot(time, v)
% 
% title('Velocity (m/s)')
% 
% 
% subplot(5,5,12)
% plot(time, M1)
% title('Mach no')
% 
% subplot(5,5,13)
% plot(time, q1)
% title('Dynamic Pressure (pa)')
% 
% subplot(5,5,14)
% hold on
% plot(time, rad2deg(gamma))
% 
% title('Trajectory Angle (Deg)')
% 
% 
% 
% subplot(5,5,15)
% plot(time, Fd1)
% title('Drag Force')
% 
% subplot(5,5,16)
% hold on
% plot(time, mFuel + 8755.1 - 994)
% title('Vehicle Mass (kg)')
% 
% 
% 
% subplot(5,5,17)
% plot(time, T1)
% title('Thrust (N)')
% 
% % Isp1 = T1./Fueldt1./9.81;
% IspNet1 = (T1-Fd1)./Fueldt1./9.81;
% 
% subplot(5,5,18)
% plot(time, Isp1)
% title('Isp')
% 
% subplot(5,5,19)
% plot(time, IspNet1)
% title('Net Isp')
% 
% subplot(5,5,20)
% plot(time, flapdeflection)
% title('Flap Deflection (deg)')
% 
% subplot(5,5,21)
% plot(time, rad2deg(Alpha))
% title('Angle of Attack (deg)')
% 
% subplot(5,5,22)
% plot(time, rad2deg(eta))
% title('Bank Angle (deg)')
% 
% subplot(5,5,23)
% plot(time, q11)
% title('Dynamic pressure after shock')
% 
% % subplot(5,5,22);
% % plot(time, dual.dynamics);
% % title('costates')
% % xlabel('time');
% % ylabel('costates');
% % legend('\lambda_1', '\lambda_2', '\lambda_3');
% 
% % subplot(5,5,23)
% % Hamiltonian = dual.Hamiltonian(1,:);
% % plot(time,Hamiltonian);
% % title('Hamiltonian')
% 
% % subplot(5,5,24)
% % hold on
% % plot(time, rad2deg(gammadot))
% % title('Trajectory Angle Change Rate (Deg/s)')
% % 
% % subplot(5,5,25)
% % hold on
% % plot(time, rad2deg(omegadot))
% % title('Omegadot Control (Deg/s2)')
% 
% 
% dim = [.8 .0 .2 .2];
% annotation('textbox',dim,'string',{['Third Stage Thrust: ', num2str(50), ' kN'],['Third Stage Starting Mass: ' num2str(2850) ' kg'],['Third Stage Isp: ' num2str(350) ' s']},'FitBoxToText','on');  
% 
% figure(202)
% sp1 = subplot(2,6,[1,6]);
% ax1 = gca; % current axes
% hold on
% plot(H/1000, alt/1000,'Color','k')
% 
% title('Trajectory')
% xlabel('Earth Normal Distance Flown (km)')
% ylabel('Vertical Position (km)')
% 
% for i = 1:floor(time(end)/30)
%     [j,k] = min(abs(time-30*i));
%     str = strcat(num2str(round(time(k))), 's');
%     text(H(k)/1000,alt(k)/1000,str,'VerticalAlignment','top', 'FontSize', 10);
%     
%     plot(H(k)/1000, alt(k)/1000, '+', 'MarkerSize', 10, 'MarkerEdgeColor','k')
% end
% 
% plot(H(end)/1000, alt(end)/1000, 'o', 'MarkerSize', 10, 'MarkerEdgeColor','k')
% 
% text(H(end)/1000,alt(end)/1000,'Third Stage Transition Point','VerticalAlignment','top', 'FontSize', 10);
% 
% dim = [.65 .45 .2 .2];
% annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass,4), ' kg'],['Second Stage Fuel Used: ' num2str(mFuel(1) - mFuel(end)) ' kg']},'FitBoxToText','on');  
% 
% thirdstageexample_H = [0+H(end) (H(end)-H(end - 1))+H(end) 20*(H(end)-H(end - 1))+H(end) 40*(H(end)-H(end - 1))+H(end) 60*(H(end)-H(end - 1))+H(end) 80*(H(end)-H(end - 1))+H(end)]/1000; %makes a small sample portion of an arbitrary third stage trajectory for example
% thirdstageexample_V = [0+alt(end) (alt(end)-alt(end - 1))+alt(end) 20*((alt(end)-alt(end -1)))+alt(end) 40*((alt(end)-alt(end -1)))+alt(end) 60*((alt(end)-alt(end -1)))+alt(end) 80*((alt(end)-alt(end -1)))+alt(end)]/1000;
% plot(thirdstageexample_H, thirdstageexample_V, 'LineStyle', '--','Color','k');
% 
% hold on
% sp2 = subplot(2,6,[7,9]);
% xlabel('time (s)')
% 
% hold on
% ax2 = gca; % current axes
% xlim([min(time) max(time)]);
% 
% line(time, rad2deg(gamma),'Parent',ax2,'Color','k', 'LineStyle','-')
% 
% line(time, M1,'Parent',ax2,'Color','k', 'LineStyle','--')
% 
% line(time, v./(10^3),'Parent',ax2,'Color','k', 'LineStyle','-.')
% 
% line(time, q1./(10^4),'Parent',ax2,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% 
% % line(time, heating_rate./(10^5),'Parent',ax1,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% % 
% % line(time, Q./(10^7),'Parent',ax1,'Color','k', 'LineStyle','-', 'lineWidth', 2.0)
% 
% % legend(ax1,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)',  'Q (Mj x 10)')
% h = legend(ax2,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)');
% rect1 = [0.12, 0.35, .25, .25];
% set(h, 'Position', rect1)
% 
% 
% sp3 = subplot(2,6,[10,12]);
% xlabel('time (s)')
% ax3 = gca;
% xlim([min(time) max(time)]);
% line(time, [rad2deg(Alpha(1:end-1)) rad2deg(Alpha(end-1))],'Parent',ax3,'Color','k', 'LineStyle','-')
% 
% line(time, flapdeflection,'Parent',ax3,'Color','k', 'LineStyle','--')
% 
% 
% % line(time, mfuel./(10^2),'Parent',ax2,'Color','k', 'LineStyle','-.')
% % line(time, eq.*10,'Parent',ax3,'Color','k', 'LineStyle','-.')
% 
% line(time, IspNet1./(10^2),'Parent',ax3,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% % 
% % g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');
% g = legend(ax3, 'AoA (degrees)','Flap Deflection (degrees)', 'Equivalence Ratio x 10', 'Net Isp (s x 10^2)');
% 
% rect2 = [0.52, 0.35, .25, .25];
% set(g, 'Position', rect2)
% 
% saveas(figure(202),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'SecondStage.fig']);
% 
% 
% 
% 
% %% PLOT RETURN
% addpath('..\SecondStageReturn\addaxis')
% 
% figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% % subplot(3,1,1)
% hold on
% 
% 
%  plot(time2,alt2/1000,'-','color','k', 'linewidth', 1.);
% % ylim([-30,40]);
% ylabel('altitude(km)');
% xlabel('time (s)');
% addaxis(time2,v2/1000,'--','color','k', 'linewidth', 1.);
% addaxislabel(2,'Velocity (km/s)');
% 
% % addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% % addaxislabel(3,'Trajectory Angle (deg)');
% 
% addaxis(time2,zeta2,':','color','k', 'linewidth', 1.2);
% addaxislabel(3,'Heading Angle (Deg)');
% 
% 
% legend(  'Altitude', 'Velocity', 'Heading Angle', 'location', 'best');
% 
% figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% % subplot(3,1,2)
% hold on
% plot(time2,rad2deg(Alpha2),'-','color','k', 'linewidth', 1.);
% ylabel('Angle of Attack (deg)');
% xlabel('time (s)')
% throttle(M2<5.0)=0;
% addaxis(time2,throttle2*100,'--','color','k', 'linewidth', 1.);
% addaxislabel(2,'Throttle (%)');
% 
% % addaxis(time,mfuel,'-.','color','k', 'linewidth', 1.);
% % addaxislabel(3,'Fuel Mass (kg)');
% 
% 
% 
% addaxis(time2,rad2deg(eta2),':','color','k', 'linewidth', 1.2);
% addaxislabel(3,'Bank Angle (Deg)');
% 
% addaxis(time2,flapdeflection2,'-.','color','k', 'linewidth', 1.2);
% addaxislabel(4,'Flap Deflection (Deg)');
% legend(  'Angle of Attack', 'Throttle' , 'Bank Angle','FlapDeflection');
% 
% 
% figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% hold on
% % subplot(3,1,3)
% plot(time2,M2,'-','color','k', 'linewidth', 1.);
% ylabel('Mach no.')
% xlabel('time (s)')
% 
% addaxis(time2,Isp2,'--','color','k', 'linewidth', 1.);
% addaxislabel(2,'Specific Impulse (s)');
% 
% addaxis(time2,q2,':','color','k', 'linewidth', 1.2);
% addaxislabel(3,'Dynamic Pressure (kPa)');
% 
% % addaxis(time,L./D,':','color','k', 'linewidth', 1.);
% % addaxislabel(4,'L/D');
% 
% legend(  'Mach no.', 'Isp (potential)', 'q' );
% 
% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % FORWARD SIMULATION
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 


forward0 = [alt(1),gamma(1),v(1),zeta(1),lat(1),lon(1)];

[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time,Alpha,f_t),ControlInterp(time,eta,f_t)),time(1:end),forward0);

figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(time,alt);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time,gamma);


subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,3));
plot(time,v);

subplot(7,1,5)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time,zeta);


subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,5));
plot(time,lat);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,6));
plot(time,lon);
% 
% % Return Forward
% forward0 = [alt2(1),gamma2(1),v2(1),zeta2(1),lat2(1),lon2(1), mFuel2(1)];
% 
% % [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
% [f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time2,Alpha2,f_t),ControlInterp(time2,eta2,f_t),ControlInterp(time2,throttle2,f_t)),time2(1:end),forward0);
% 
% % altitude  = (output.result.solution.phase(1).state(:,1)-auxdata.Re);
% figure(213)
% subplot(7,1,1)
% hold on
% plot(f_t(1:end),f_y(:,1));
% plot(time2,alt2);
% 
% % gamma  = output.result.solution.phase.state(:,5);
% 
% subplot(7,1,2)
% hold on
% plot(f_t(1:end),f_y(:,2));
% plot(time2,gamma2);
% 
% % latitude  = output.result.solution.phase.state(:,3);
% subplot(7,1,3:5)
% hold on
% plot(f_y(:,6),f_y(:,5));
% plot(lon2,lat2);
% 
% subplot(7,1,6)
% hold on
% plot(f_t(1:end),f_y(:,7));
% plot(time2,mFuel2);
% 
% %% Check KKT and pontryagins minimum
% % Check that the hamiltonian = 0 (for free end time)
% % Necessary condition
% input_test = output.result.solution;
% input_test.auxdata = auxdata;
% phaseout_test = CombinedContinuous(input_test);
% 
% lambda1 = output.result.solution.phase(1).costate;
% for i = 1:length(lambda1)-1
%     H1(i) = lambda1(i+1,:)*phaseout_test(1).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
% end
% 
% lambda2 = output.result.solution.phase(2).costate;
% for i = 1:length(lambda2)-1
%     H2(i) = lambda2(i+1,:)*phaseout_test(2).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
% end
% 
% figure(221)
% hold on
% plot(time(1:end-1),H1)
% plot(time2(1:end-1),H2)
% ylabel('Hamiltonian')
% xlabel('Time (s)')
% legend('Ascent','Return')
% 
% % Check Primal Feasibility
% % Check calculated derivatives with the numerical derivative of each
% % porimal, scaled by that primal
% figure(220)
% hold on
% for i = 1:length(output.result.solution.phase(1).state(1,:))
% plot(time,([diff(output.result.solution.phase(1).state(:,i))./diff(output.result.solution.phase(1).time); 0] - phaseout_test(1).dynamics(:,i))./output.result.solution.phase(1).state(:,i),'--');
% end
% for i = 1:length(output.result.solution.phase(2).state(1,:))
%     if i<= 7 % Plot different line styles when no. of colours exceeded
%     plot(time2,([diff(output.result.solution.phase(2).state(:,i))./diff(output.result.solution.phase(2).time); 0] - phaseout_test(2).dynamics(:,i))./output.result.solution.phase(2).state(:,i));
%     else
%     plot(time2,([diff(output.result.solution.phase(2).state(:,i))./diff(output.result.solution.phase(2).time); 0] - phaseout_test(2).dynamics(:,i))./output.result.solution.phase(2).state(:,i),':');
%     end
% end
% xlabel('Time (s)')
% ylabel('Derivative Error')
% ylim([-1,1])
% legend('Alt Ascent','lon Ascent','lat Ascent','v Ascent','gamma Ascent','zeta Ascent','aoa Ascent','bank Ascent','mFuel Ascent', 'Alt Descent','lon Descent','lat Descent','v Descent','gamma Descent','zeta Descent','aoa Descent','bank Descent','mFuel Descent','throttle Descent')
% 
% %% plot engine interpolation visualiser
% T0 = spline( auxdata.interp.Atmosphere(:,1),  auxdata.interp.Atmosphere(:,2), alt); 
% T_in1 = auxdata.interp.tempgridded(M1,rad2deg(Alpha)).*T0;
% M_in1 = auxdata.interp.M1gridded(M1, rad2deg(Alpha));
% 
% plotM = [min(M_englist):0.01:9];
% plotT = [min(T_englist):1:550];
% [gridM,gridT] =  ndgrid(plotM,plotT);
% interpeq = auxdata.interp.eqGridded(gridM,gridT);
% interpIsp = auxdata.interp.IspGridded(gridM,gridT);
% 
% figure(210)
% hold on
% contourf(gridM,gridT,interpeq,50);
% scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,4),'filled');
% xlabel('M1')
% ylabel('T1')
% plot(M_in1,T_in1,'r');
% 
% error_Isp = auxdata.interp.IspGridded(engine_data(:,1),engine_data(:,2))-engine_data(:,3);
% 
% figure(211)
% hold on
% contourf(gridM,gridT,interpIsp,100,'LineWidth',0);
% scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,3),'k')
% xlabel('M1')
% ylabel('T1')
% c=colorbar
% c.Label.String = 'ISP';
% plot(M_in1,T_in1,'r');
% 
% %%
% [gridM2,gridAoA2] =  ndgrid(plotM,plotT);
% 
% 
% 
% % Run First Stage =========================================================
% const_firststage = 1;
% addpath('../../FirstStage')
% % addpath('../../DIDO_7.3.7')
% % run startup.m
% [FirstStageStates] = FirstStageProblem(alt(1),gamma(1),lat(1),zeta(1),const_firststage);
% cd('../SecondStage/Combined - 2ns Stage Ascent and Return')
% dlmwrite('FirstStage.txt', FirstStageStates);
% copyfile('FirstStage.txt',sprintf('../ArchivedResults/%s/firststage_%s.txt',Timestamp,Timestamp))
% 
% 
% %% Latitude Plot
% figure(250)
% plot(FirstStageStates(:,9))
% plot(phi)
% plot(ThirdStagePhi)
% title('Latitude')
% 
% %% SAVE FIGS
% saveas(figure(301),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ThirdStage.fig']);
% saveas(figure(101),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'FirstStage.fig']);
% %%
% 
% % =========================================================================
% % Troubleshooting Procedure
% % =========================================================================
% 
% % 1: Check that you have posed your problem correctly ie. it is physically
% % feasible and the bounds allow for a solution
% % 2: Check for NaN values (check derivatives in Dynamics file while running)
% % 3: Check guess, is it reasonable? Is it too close to the expected
% % solution? Both can cause errors! Sometimes there is no real rhyme or
% % reason to picking the correct guess, but a close bound to
% % the expected solution has worked the most in my experience
% % 4: Play with the no. of nodes, try both even and odd values
% % 5: Play with scaling
% % 6: Try all of the above in various combinations until it works!
% 
% 



