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

Stage1.A = 93.25; %Reference Area in m² 
Stage1.mStruct = 205245.8;

Stage1.T_SL = 1961*9*1e3; %for mixture ratio 6 with all boosters active
Stage1.Isp_SL = 389;

Stage1.T_vac = 2206*9*1e3;
Stage1.Isp_vac = 437;

Stage2.A = 461; %Reference Area in m² 
Stage2.mStruct = 134361.1;

Stage2.T_SL = 1830*2*1e3; %for mixture ratio 6 with all boosters active
Stage2.Isp_SL = 363;

Stage2.T_vac = 2268*2*1e3;
Stage2.Isp_vac = 449;

auxdata.Stage1 = Stage1;
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
Stage1_MList = [0.00  0.20  0.40  0.60  0.80  1.00  1.20  1.40  1.60  1.80  2.00  2.20  2.40  2.60  2.80  3.00  3.20  3.40  3.60  3.80  4.00 5. 6. 7. 8. 9. 10. 11. 12. 13. 14. 15. 18. 21. 25.];
Stage1_aoaList = [ 0.00   1.00   2.00   3.00   4.00   5.00   6.00   7.00   8.00   9.00  10.00  11.00  12.00  13.00  14.00  15.00  16.00  17.00  18.00  19.00  20.00];
Stage1_Cl = importdata('aero_BoosterandOrbiter_Lift');
Stage1_Cd = importdata('aero_BoosterandOrbiter_Drag');

Stage1_MList_Full = zeros(length(Stage1_Cl),1);
Stage1_aoaList_Full = zeros(length(Stage1_Cl),1);

for i = 1:length(Stage1_Cl)
    Stage1_aoaList_Full(i) = Stage1_aoaList(i - length(Stage1_aoaList)*floor((i-1)/length(Stage1_aoaList))); 
    Stage1_MList_Full(i) = Stage1_MList(floor((i-1)/length(Stage1_aoaList)) + 1);
end

[Stage1_MGrid,Stage1_aoaGrid] = meshgrid(Stage1_MList,Stage1_aoaList);

Stage1_Cl_scattered = scatteredInterpolant(Stage1_MList_Full,Stage1_aoaList_Full,Stage1_Cl);
Stage1_Cd_scattered = scatteredInterpolant(Stage1_MList_Full,Stage1_aoaList_Full,Stage1_Cd);

for i = 1:length(Stage1_MList)
    for j = 1:length(Stage1_aoaList)
        Stage1_Cl_Grid(j,i) = Stage1_Cl_scattered(Stage1_MList(i),Stage1_aoaList(j));
        Stage1_Cd_Grid(j,i) = Stage1_Cd_scattered(Stage1_MList(i),Stage1_aoaList(j));
    end
end

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
        Stage2_Cl_Grid(j,i) = Stage2_Cl_scattered(Stage2_MList(i),Stage2_aoaList(j));
        Stage2_Cd_Grid(j,i) = Stage2_Cd_scattered(Stage2_MList(i),Stage2_aoaList(j));
    end
end

%% Create Gridded Interpolants

auxdata.interp.Stage1.Cl_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cl_Grid','linear', 'linear');
auxdata.interp.Stage1.Cd_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cd_Grid','linear', 'linear');

auxdata.interp.Stage2.Cl_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cl_Grid','linear', 'linear');
auxdata.interp.Stage2.Cd_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cd_Grid','linear', 'linear');

%% Get population interpolation
load PopInterp

auxdata.PopInterp = PopInterp;

%% Import Bounds %%========================================================
lonMin = -2*pi;         
lonMax = 2*pi;

latMin = -pi/2+0.0000001;  
% latMin = 0; 
latMax = pi/2-0.0000001;


%% Initial Geo

lat0 = deg2rad(-23.3791); % Rockhampton
lon0 = deg2rad(150.5100); % Rockhampton

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

% lat0 = deg2rad(45.37); %South Japan after true launch
% lon0 = deg2rad(137.71);

auxdata.lon0 = lon0;
%%
aoaMin = 0;  
aoaMax = deg2rad(25);

bankMin = deg2rad(-50); 
bankMax =   deg2rad(50);

% aoaMin = 0;  aoaMax = 20*pi/180;
% bankMin1 = -50*pi/180; bankMax1 =   50*pi/180;

% Primal Bounds
bounds.phase(1).state.lower = [0, lonMin, latMin, 0, -deg2rad(50), -pi, 0, aoaMin, bankMin];
bounds.phase(1).state.upper = [150000, lonMax, latMax, 15000, deg2rad(85), pi, 1.6038e+06, aoaMax, bankMax];

% Initial States
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 95, deg2rad(80), deg2rad(70), 1.5038e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 105, deg2rad(89), deg2rad(80), 1.6038e+06, aoaMax];

bounds.phase(1).initialstate.lower = [1000,0, lat0, 100, deg2rad(80), 0, 1.28e+06, aoaMin, bankMin] ;
bounds.phase(1).initialstate.upper = [1200,0, lat0, 120, deg2rad(85), 0, 1.32e+06, aoaMax, bankMax];

bounds.phase(1).finalstate.lower = bounds.phase(1).state.lower;
bounds.phase(1).finalstate.upper = bounds.phase(1).state.upper;

bounds.phase(2).initialstate.lower = bounds.phase(1).state.lower;
bounds.phase(2).initialstate.upper = bounds.phase(1).state.upper;

% Control Bounds
% bounds.phase(1).control.lower = [deg2rad(-.5)];
% bounds.phase(1).control.upper = [deg2rad(.5)];

bounds.phase(1).control.lower = [deg2rad(-.5), deg2rad(-1)];
bounds.phase(1).control.upper = [deg2rad(.5), deg2rad(1)];


% Time Bounds
bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower = 50;
bounds.phase(1).finaltime.upper = 5000;

%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.

bounds.phase(1).path.lower = [0, -2.5*9.81];
bounds.phase(1).path.upper = [60000, 2.5*9.81];

%% Bound integral if necessary
bounds.phase(1).integral.lower = 0;
bounds.phase(1).integral.upper = 1e50;


%% Set all phase bounds
bounds.phase(2).state = bounds.phase(1).state;
bounds.phase(2).finalstate = bounds.phase(1).finalstate;

bounds.phase(2).path = bounds.phase(1).path;

bounds.phase(2).integral = bounds.phase(1).integral;

bounds.phase(2).initialtime.lower = 0;
bounds.phase(2).initialtime.upper = 5000;

bounds.phase(2).finaltime = bounds.phase(1).finaltime;

bounds.phase(2).control = bounds.phase(1).control;

bounds.phase(3) = bounds.phase(2);
bounds.phase(4) = bounds.phase(2);
bounds.phase(5) = bounds.phase(2);
bounds.phase(6) = bounds.phase(2);
bounds.phase(7) = bounds.phase(2);
bounds.phase(8) = bounds.phase(2);

bounds.phase(8).finalstate.lower =  [70000, lonMin, latMin, 7000, 0, -pi, 0, aoaMin, bankMin];
bounds.phase(8).finalstate.upper = [70000, lonMax, latMax, 7000, 0, pi, 1.6038e+06, aoaMax, bankMax];



% bounds.phase(8).initialtime = bounds.phase(2).initialtime;
% bounds.phase(8).finaltime = bounds.phase(2).finaltime;
% bounds.phase(8).state = bounds.phase(2).state;
% bounds.phase(8).initialstate = bounds.phase(2).initialstate;
% bounds.phase(8).control = bounds.phase(2).control;
% End States

% bounds.phase(8).finalstate.lower = [70000, lonMin, latMin, 0, 0, -2*pi, 0, aoaMin];
% bounds.phase(8).finalstate.upper = [150000, lonMax, latMax, 15000, deg2rad(80), 2*pi,1.6038e+06, aoaMax];
% 
% bounds.phase(8).finalstate.lower = [70000, lonMin, latMin, 0, 0, -2*pi, 0, aoaMin];
% bounds.phase(8).finalstate.upper = [70000, lonMax, latMax, 15000, deg2rad(0), 2*pi,1.6038e+06, aoaMax];


%% Event bounds
bounds.eventgroup(1).lower = zeros(1,10);
bounds.eventgroup(1).upper = zeros(1,10);

bounds.eventgroup(2).lower = zeros(1,10);
bounds.eventgroup(2).upper = zeros(1,10);

bounds.eventgroup(3).lower = zeros(1,10);
bounds.eventgroup(3).upper = zeros(1,10);

bounds.eventgroup(4).lower = zeros(1,10);
bounds.eventgroup(4).upper = zeros(1,10);

bounds.eventgroup(5).lower = zeros(1,10);
bounds.eventgroup(5).upper = zeros(1,10);

bounds.eventgroup(6).lower = zeros(1,10);
bounds.eventgroup(6).upper = zeros(1,10);

bounds.eventgroup(7).lower = zeros(1,10);
bounds.eventgroup(7).upper = zeros(1,10);

bounds.eventgroup(8).lower = zeros(1,10);
bounds.eventgroup(8).upper = zeros(1,10);

bounds.eventgroup(9).lower = ones(1,9);
bounds.eventgroup(9).upper = [500*ones(1,8) 10000];

%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(1).state(:,1)   = [2000;7000];
% guess.phase(1).state(:,2)   = [2.5;2.55];
guess.phase(1).state(:,2)   = [0;0.05];
guess.phase(1).state(:,3)   = [lat0;lat0+0.05];
guess.phase(1).state(:,4)   = [100,1100];
guess.phase(1).state(:,5)   = [deg2rad(80),deg2rad(80)];
guess.phase(1).state(:,6)   = [deg2rad(0),deg2rad(0)];
guess.phase(1).state(:,7) 	= [1.6038e+06, 1.500e+06];
guess.phase(1).state(:,8)   = [1*pi/180; 1*pi/180];
guess.phase(1).state(:,9)   = [deg2rad(0);deg2rad(0)];


guess.phase(1).control     = [[0;0],[0;0]];
guess.phase(1).time          = [0;50];


guess.phase(1).integral = 0;
% Tire stages together
% bounds.eventgroup(1).lower = [zeros(1,10)];
% bounds.eventgroup(1).upper = [zeros(1,10)]; 

guess.phase(2).state = guess.phase(1).state + [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0;0 0]';
guess.phase(3).state = guess.phase(2).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0;0 0]';
guess.phase(4).state = guess.phase(3).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0;0 0]';
guess.phase(5).state = guess.phase(4).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0;0 0]';
guess.phase(6).state = guess.phase(5).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0;0 0]';
guess.phase(7).state = guess.phase(6).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0;0 0]';
guess.phase(8).state = guess.phase(7).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0;0 0]';

guess.phase(2).control = guess.phase(1).control;
guess.phase(3).control = guess.phase(2).control;
guess.phase(4).control = guess.phase(3).control;
guess.phase(5).control = guess.phase(4).control;
guess.phase(6).control = guess.phase(5).control;
guess.phase(7).control = guess.phase(6).control;
guess.phase(8).control = guess.phase(7).control;

guess.phase(2).time = guess.phase(1).time + 50;
guess.phase(3).time = guess.phase(2).time + 50;
guess.phase(4).time = guess.phase(3).time + 50;
guess.phase(5).time = guess.phase(4).time + 50;
guess.phase(6).time = guess.phase(5).time + 50;
guess.phase(7).time = guess.phase(6).time + 50;
guess.phase(8).time = guess.phase(7).time + 50;

guess.phase(2).integral = guess.phase(1).integral;
guess.phase(3).integral = guess.phase(2).integral;
guess.phase(4).integral = guess.phase(3).integral;
guess.phase(5).integral = guess.phase(4).integral;
guess.phase(6).integral = guess.phase(5).integral;
guess.phase(7).integral = guess.phase(6).integral;
guess.phase(8).integral = guess.phase(7).integral;



%% DESCENT BOUNDS

altMin = 1;
altMax = 150000;


vMin = 1;
vMax = 15000;

gammaMin = -deg2rad(80);
gammaMax = deg2rad(80);

zetaMin = -2*pi;
% zetaMin = 0;
zetaMax = 2*pi;


% latMax = pi;


% bankMin = deg2rad(-60); 
% bankMax =   deg2rad(60);


% alt0 = 70000;
% v0 = 7000;
% gamma0 = 0;

%japan-germany
alt0 = 73235;
v0 = 6523;
gamma0 = deg2rad(0.128);

%% End Geo
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
%%

% Primal Bounds
bounds.phase(9).state.lower = [altMin, lonMin, latMin, vMin, gammaMin, zetaMin, 0, aoaMin, bankMin];
bounds.phase(9).state.upper = [altMax, lonMax, latMax, vMax, gammaMax, zetaMax, 1.6038e+06, aoaMax, bankMax];

% Initial States
% bounds.phase(1).initialstate.lower = [alt0,0, lat0, v0, gamma0, zetaMin, aoaMin, bankMin] ;
% bounds.phase(1).initialstate.upper = [alt0,0, lat0, v0, gamma0, zetaMax, aoaMax, bankMax];

% bounds.phase(1).initialstate.lower = [alt0,0, lat0, v0-1000, gamma0, zetaMin, aoaMin, bankMin] ;
% bounds.phase(1).initialstate.upper = [alt0,0, lat0, v0, gamma0, zetaMax, aoaMax, bankMax];

bounds.phase(9).initialstate.lower = [70000,lonMin, latMin, 7000, 0, zetaMin, 0, aoaMin, bankMin] ;
bounds.phase(9).initialstate.upper = [70000,lonMax, latMax, 7000, 0, zetaMax, 1.6038e+06, aoaMax, bankMax];


% End States

%Starting East heading east
% For aus-germany, japan-germany
bounds.phase(9).finalstate.lower = [altFMin, lonF-lon0+2*pi, latF, vMin, deg2rad(-20), zetaMin, 0, aoaMin, 0];
bounds.phase(9).finalstate.upper = [altFMax, lonF-lon0+2*pi, latF, 150, deg2rad(20), zetaMax, 0, aoaMax, 0];

% bounds.phase(9).finalstate.lower = [altMin, lonF-lon0+2*pi, latF, vMin, gammaMin, zetaMin, aoaMin, bankMin] ;
% bounds.phase(9).finalstate.upper = [altMax, lonF-lon0+2*pi, latF, vMax, gammaMax, zetaMax, aoaMax, bankMax];

%  bounds.phase(9).finalstate.lower = [altMin, lonMin, latMin, vMin, gammaMin, zetaMin, 0, aoaMin, bankMin] ;
% bounds.phase(9).finalstate.upper = [altMax, lonMax, latMax, vMax, gammaMax, zetaMax, 1.6038e+06, aoaMax, bankMax];

% bounds.phase(9).finalstate.lower = [altMin, lonMin, latMin, vMin, -0.01, zetaMin, 0, aoaMin, bankMin] ;
% bounds.phase(9).finalstate.upper = [3000, lonMax, latMax, 300, 0.01, zetaMax, 1.6038e+06, aoaMax, bankMax];


%Starting West heading west
% for germany-japan
% bounds.phase(1).finalstate.lower = [altF, lonF-lon0-2*pi, latF, vMin, -0.05, zetaMin, aoaMin, 0];
% bounds.phase(1).finalstate.upper = [altF, lonF-lon0-2*pi, latF, 50, 0.05, zetaMax, aoaMax, 0];

% If not going outside of degree bounds
% for japan-florida, florida-japan
% bounds.phase(1).finalstate.lower = [altF, lonF-lon0, latF, vMin, -0.05, zetaMin, aoaMin, 0];
% bounds.phase(1).finalstate.upper = [altF, lonF-lon0, latF, 50, 0.05, zetaMax, aoaMax, 0];

% Control Bounds
bounds.phase(9).control.lower = [deg2rad(-.1), deg2rad(-1)];
bounds.phase(9).control.upper = [deg2rad(.1), deg2rad(1)];

% Time Bounds
bounds.phase(9).initialtime.lower = 0;
bounds.phase(9).initialtime.upper = 10000;
bounds.phase(9).finaltime.lower =0;
bounds.phase(9).finaltime.upper = 10000;

%Integral Bounds

bounds.phase(9).integral.lower = 0;
bounds.phase(9).integral.upper = 1e50;

% path bounds
bounds.phase(1).path.lower = [0, -2.5*9.81];
bounds.phase(1).path.upper = [60000, 2.5*9.81];

%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(9).state(:,1)   = [guess.phase(8).state(2,1);50000];
guess.phase(9).state(:,2)   = [guess.phase(8).state(2,2);lonF-lon0+2*pi]; 
% guess.phase(1).state(:,2)   = [0;lonF-lon0-2*pi];
% guess.phase(1).state(:,2)   = [0;lonF-lon0];
guess.phase(9).state(:,3)   = [guess.phase(8).state(2,3);latF];
% guess.phase(1).state(:,3)   = [lat0;latF-0.5];
guess.phase(9).state(:,4)   = [guess.phase(8).state(2,4);vMin];

% This is gamma
guess.phase(9).state(:,5)   = [guess.phase(8).state(2,5);0]; %

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
% guess.phase(9).state(:,6)   = [deg2rad(80);deg2rad(-60)]; % japan-Germany actual

guess.phase(9).state(:,6)   = [deg2rad(0);deg2rad(0)];

guess.phase(9).state(:,7)   = [deg2rad(0);deg2rad(0)];


 guess.phase(9).state(:,8)   = [1*pi/180; 1*pi/180];
% guess.phase(9).state(:,8)   = [0;0];
guess.phase(9).state(:,9)   = [deg2rad(50);deg2rad(50)];

guess.phase(9).control      = [[0;0],[0;0]];
guess.phase(9).time          = [400;4000];

guess.phase(9).integral = 0;



%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
% mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 3;
mesh.colpointsmin = 3;
mesh.colpointsmax = 50;
mesh.tolerance    = 1e-5;


%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                           = 'SpaceLinerCombined';
setup.functions.continuous           = @SpaceLinerContinuousCombined;
setup.functions.endpoint             = @SpaceLinerEndpointCombined;
setup.auxdata                        = auxdata;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.displaylevel                   = 2;
setup.nlp.solver                     = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.maxiterations = 1000;
setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.derivativelevel    = 'first';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
setup.scales.method                  = 'automatic-guess';

%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%


output = gpops2(setup);

%%

EndTime = datestr(now,30) % Display the ending time

% =========================================================================

alt = [output.result.solution.phase(1).state(:,1).' ...
    output.result.solution.phase(2).state(:,1).' ...
    output.result.solution.phase(3).state(:,1).' ...
    output.result.solution.phase(4).state(:,1).' ...
    output.result.solution.phase(5).state(:,1).' ...
    output.result.solution.phase(6).state(:,1).' ...
    output.result.solution.phase(7).state(:,1).' ...
    output.result.solution.phase(8).state(:,1).' ...
    output.result.solution.phase(9).state(:,1).'];
lon = [output.result.solution.phase(1).state(:,2).' ...
    output.result.solution.phase(2).state(:,2).' ...
    output.result.solution.phase(3).state(:,2).' ...
    output.result.solution.phase(4).state(:,2).' ...
    output.result.solution.phase(5).state(:,2).' ...
    output.result.solution.phase(6).state(:,2).' ...
    output.result.solution.phase(7).state(:,2).' ...
    output.result.solution.phase(8).state(:,2).' ...
    output.result.solution.phase(9).state(:,2).'];
lat = [output.result.solution.phase(1).state(:,3).' ...
    output.result.solution.phase(2).state(:,3).' ...
    output.result.solution.phase(3).state(:,3).' ...
    output.result.solution.phase(4).state(:,3).' ...
    output.result.solution.phase(5).state(:,3).' ...
    output.result.solution.phase(6).state(:,3).' ...
    output.result.solution.phase(7).state(:,3).' ...
    output.result.solution.phase(8).state(:,3).' ...
    output.result.solution.phase(9).state(:,3).'];
v = [output.result.solution.phase(1).state(:,4).' ...
    output.result.solution.phase(2).state(:,4).' ...
    output.result.solution.phase(3).state(:,4).' ...
    output.result.solution.phase(4).state(:,4).' ...
    output.result.solution.phase(5).state(:,4).' ...
    output.result.solution.phase(6).state(:,4).' ...
    output.result.solution.phase(7).state(:,4).' ...
    output.result.solution.phase(8).state(:,4).' ...
    output.result.solution.phase(9).state(:,4).'];
gamma = [output.result.solution.phase(1).state(:,5).' ...
    output.result.solution.phase(2).state(:,5).' ...
    output.result.solution.phase(3).state(:,5).' ...
    output.result.solution.phase(4).state(:,5).' ...
    output.result.solution.phase(5).state(:,5).' ...
    output.result.solution.phase(6).state(:,5).' ...
    output.result.solution.phase(7).state(:,5).' ...
    output.result.solution.phase(8).state(:,5).' ...
    output.result.solution.phase(9).state(:,5).'];
zeta = [output.result.solution.phase(1).state(:,6).' ...
    output.result.solution.phase(2).state(:,6).' ...
    output.result.solution.phase(3).state(:,6).' ...
    output.result.solution.phase(4).state(:,6).' ...
    output.result.solution.phase(5).state(:,6).' ...
    output.result.solution.phase(6).state(:,6).' ...
    output.result.solution.phase(7).state(:,6).' ...
    output.result.solution.phase(8).state(:,6).' ...
    output.result.solution.phase(9).state(:,6).'];
mFuel = [output.result.solution.phase(1).state(:,7).' ...
    output.result.solution.phase(2).state(:,7).' ...
    output.result.solution.phase(3).state(:,7).' ...
    output.result.solution.phase(4).state(:,7).' ...
    output.result.solution.phase(5).state(:,7).' ...
    output.result.solution.phase(6).state(:,7).' ...
    output.result.solution.phase(7).state(:,7).' ...
    output.result.solution.phase(8).state(:,7).' ...
    output.result.solution.phase(9).state(:,7).']; 
Alpha = [output.result.solution.phase(1).state(:,8).' ...
    output.result.solution.phase(2).state(:,8).' ...
    output.result.solution.phase(3).state(:,8).' ...
    output.result.solution.phase(4).state(:,8).' ...
    output.result.solution.phase(5).state(:,8).' ...
    output.result.solution.phase(6).state(:,8).' ...
    output.result.solution.phase(7).state(:,8).' ...
    output.result.solution.phase(8).state(:,8).' ...
    output.result.solution.phase(9).state(:,8).'];

eta = [output.result.solution.phase(1).state(:,9).' ...
    output.result.solution.phase(2).state(:,9).' ...
    output.result.solution.phase(3).state(:,9).' ...
    output.result.solution.phase(4).state(:,9).' ...
    output.result.solution.phase(5).state(:,9).' ...
    output.result.solution.phase(6).state(:,9).' ...
    output.result.solution.phase(7).state(:,9).' ...
    output.result.solution.phase(8).state(:,9).' ...
    output.result.solution.phase(9).state(:,9).'];



% eta = output.result.solution.phase(1).state(:,9).';

states.alt     = output.result.solution.phase(1).state(:,1);
states.lon     = output.result.solution.phase(1).state(:,2);
states.lat     = output.result.solution.phase(1).state(:,3);
states.v       = output.result.solution.phase(1).state(:,4);
states.gamma   = output.result.solution.phase(1).state(:,5);
states.zeta    = output.result.solution.phase(1).state(:,6);
states.mFuel   = output.result.solution.phase(1).state(:,7);



time = [output.result.solution.phase(1).time.' ...
    output.result.solution.phase(2).time.' ...
    output.result.solution.phase(3).time.' ...
    output.result.solution.phase(4).time.' ...
    output.result.solution.phase(5).time.' ...
    output.result.solution.phase(6).time.' ...
    output.result.solution.phase(7).time.' ...
    output.result.solution.phase(8).time.' ...
    output.result.solution.phase(9).time.'];

time_ascent = [output.result.solution.phase(1).time.' ...
    output.result.solution.phase(2).time.' ...
    output.result.solution.phase(3).time.' ...
    output.result.solution.phase(4).time.' ...
    output.result.solution.phase(5).time.' ...
    output.result.solution.phase(6).time.' ...
    output.result.solution.phase(7).time.' ...
    output.result.solution.phase(8).time.'];

time_descent = output.result.solution.phase(9).time.';

% eta = [zeros(1,length(time_ascent)) output.result.solution.phase(9).state(:,8).' ];



controls.Alpha = Alpha';




% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 

throttle = ones(length(alt),1);

throttle(time>output.result.solution.phase(1).time(end)) = 0.891;

throttle(time>output.result.solution.phase(2).time(end)) = 0.812;

throttle(time>output.result.solution.phase(3).time(end)) = .7333;

throttle(time>output.result.solution.phase(4).time(end)) = .6545;

throttle(time>output.result.solution.phase(5).time(end)) = .5757;

throttle(time>output.result.solution.phase(6).time(end)) = 0.496;

throttle(time>=output.result.solution.phase(7).time(end)) = 1; %after separation

throttle(time>=output.result.solution.phase(8).time(end)) = 0; %after separation

 
auxdata.mFuel_descent = output.result.solution.phase(8).state(end,7);

forward0 = [alt(1),gamma(1),v(1),zeta(1),lat(1),lon(1), mFuel(1)];


time_diff = [1 diff(time)];

time_forward = time;
Alpha_forward = Alpha;
eta_forward = eta;
throttle_forward = throttle;

time_forward(time_diff==0) = [];
Alpha_forward(time_diff==0) = [];
eta_forward(time_diff==0) = [];
throttle_forward(time_diff==0) = [];

stage_forward = ones(1,length(time_forward));
stage_forward(time_forward>=output.result.solution.phase(7).time(end)) = 2;
% stage_forward(time_forward>=output.result.solution.phase(8).time(end)) = 3;


[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time_forward,Alpha_forward,f_t),0,ControlInterp(time_forward,throttle_forward,f_t),ControlInterp(time_forward,stage_forward,f_t)),0:time_ascent(end),forward0);


forward0_descent = f_y(end,1:6);
[f_t2, f_y2] = ode45(@(f_t2,f_y2) VehicleModel_forward(f_t2, f_y2,auxdata,ControlInterp(time_forward,Alpha_forward,f_t2),ControlInterp(time_forward,eta_forward,f_t2),0,3),time_descent(2):time_descent(end),forward0_descent);



% [f_t, f_y] = ode23(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time,Alpha,f_t),ControlInterp(time,throttle,f_t),time(end)),time(1:end),forward0);

% 
% dt = 1;
% f_y = forward0;
% temp = 1;
% for t_temp = 0:dt:time(end)
%     
%     df_y = VehicleModel_forward(t_temp, f_y(temp,:),auxdata,ControlInterp(time,Alpha,t_temp),ControlInterp(time,eta,t_temp),ControlInterp(time,throttle,t_temp),ControlInterp(time_forward,stage_forward,f_t));
%   
%     
%     f_y(temp+1,:) = f_y(temp,:) + dt*df_y';  
%     temp = temp+1;
% end
% f_t = 0:dt:time(end);


figure(212)
subplot(7,1,[1 2])
hold on
plot(f_t(1:end),f_y(:,1));
plot(f_t2(1:end),f_y2(:,1));
plot(time,alt);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(f_t2(1:end),f_y2(:,2));
plot(time,gamma);


subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,3));
plot(f_t2(1:end),f_y2(:,3));
plot(time,v);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,4));
plot(f_t2(1:end),f_y2(:,4));
plot(time,zeta);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,7));
% plot(f_t2(1:end),f_y2(:,7));
plot(time,mFuel);

figure()
hold on
plot(f_y(:,6),f_y(:,5));
plot(f_y2(:,6),f_y2(:,5));
plot(lon,lat);


phase.state(:,1) = alt;
phase.state(:,2) = lon;
phase.state(:,3) = lat;
phase.state(:,4) = v;
phase.state(:,5) = gamma;
phase.state(:,6) = zeta;
phase.state(:,7) = mFuel;
phase.state(:,8) = Alpha;
phase.state(:,9) = eta;

stage = ones(1,length(time));
stage(time>=output.result.solution.phase(7).time(end)) = 2;

T = [];
M = [];
Fd = [];
L = [];
q1 = [];
m = [];
for i = 1:length(alt)
    phase_temp.state = phase.state(i,:);
[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1(i), M(i), Fd(i), rho,L(i),Fueldt1,T(i),Isp1,Isp2,m(i)] = SpaceLinerVehicleModel(time(i),phase_temp,throttle(i),auxdata,stage(i));

end

figure(201)
subplot(5,2,1)
hold on
plot(time,alt)
xlabel('time')
ylabel('altitude')

subplot(5,2,2)
hold on
plot(v,alt)
xlabel('velocity')
ylabel('altitude')


subplot(5,2,3)
hold on
plot(time,rad2deg(Alpha))
plot(time,rad2deg(gamma))

subplot(5,2,4)
hold on
plot(time,T)
plot(time,Fd)
plot(time,L)
xlabel('time')

subplot(5,2,6)
hold on
plot(time,m)

subplot(5,2,7)
hold on
plot(time,M)

subplot(5,2,8)
hold on
plot(time,q1)

figure(230)
hold on
plot3(lon,lat,alt)


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

%% Check KKT and pontryagins minimum
% Check that the hamiltonian = 0 (for free end time)
% Necessary condition
input_test = output.result.solution;
input_test.auxdata = auxdata;
phaseout_test = CombinedContinuous(input_test);

lambda1 = output.result.solution.phase(1).costate;
for i = 1:length(lambda1)-1
    H1(i) = lambda1(i+1,:)*phaseout_test(1).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

lambda2 = output.result.solution.phase(2).costate;
for i = 1:length(lambda2)-1
    H2(i) = lambda2(i+1,:)*phaseout_test(2).dynamics(i,:).'; %H = lambda transpose * f(x,u,t) + L, note that there is no continuous cost L
end

figure(221)
hold on
plot(time(1:end-1),H1)
plot(time2(1:end-1),H2)
ylabel('Hamiltonian')
xlabel('Time (s)')
legend('Ascent','Return')

% Check Primal Feasibility
% Check calculated derivatives with the numerical derivative of each
% porimal, scaled by that primal
figure(220)
hold on
for i = 1:length(output.result.solution.phase(1).state(1,:))
plot(time,([diff(output.result.solution.phase(1).state(:,i))./diff(output.result.solution.phase(1).time); 0] - phaseout_test(1).dynamics(:,i))./output.result.solution.phase(1).state(:,i),'--');
end
for i = 1:length(output.result.solution.phase(2).state(1,:))
    if i<= 7 % Plot different line styles when no. of colours exceeded
    plot(time2,([diff(output.result.solution.phase(2).state(:,i))./diff(output.result.solution.phase(2).time); 0] - phaseout_test(2).dynamics(:,i))./output.result.solution.phase(2).state(:,i));
    else
    plot(time2,([diff(output.result.solution.phase(2).state(:,i))./diff(output.result.solution.phase(2).time); 0] - phaseout_test(2).dynamics(:,i))./output.result.solution.phase(2).state(:,i),':');
    end
end
xlabel('Time (s)')
ylabel('Derivative Error')
ylim([-1,1])
legend('Alt Ascent','lon Ascent','lat Ascent','v Ascent','gamma Ascent','zeta Ascent','aoa Ascent','bank Ascent','mFuel Ascent', 'Alt Descent','lon Descent','lat Descent','v Descent','gamma Descent','zeta Descent','aoa Descent','bank Descent','mFuel Descent','throttle Descent')

%% plot engine interpolation visualiser
T0 = spline( auxdata.interp.Atmosphere(:,1),  auxdata.interp.Atmosphere(:,2), alt); 
T_in1 = auxdata.interp.tempgridded(M1,rad2deg(Alpha)).*T0;
M_in1 = auxdata.interp.M1gridded(M1, rad2deg(Alpha));

plotM = [min(M_englist):0.01:9];
plotT = [min(T_englist):1:550];
[gridM,gridT] =  ndgrid(plotM,plotT);
interpeq = auxdata.interp.eqGridded(gridM,gridT);
interpIsp = auxdata.interp.IspGridded(gridM,gridT);

figure(210)
hold on
contourf(gridM,gridT,interpeq,50);
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,4),'filled');
xlabel('M1')
ylabel('T1')
plot(M_in1,T_in1,'r');

error_Isp = auxdata.interp.IspGridded(engine_data(:,1),engine_data(:,2))-engine_data(:,3);

figure(211)
hold on
contourf(gridM,gridT,interpIsp,100,'LineWidth',0);
scatter(engine_data(:,1),engine_data(:,2),30,engine_data(:,3),'k')
xlabel('M1')
ylabel('T1')
c=colorbar
c.Label.String = 'ISP';
plot(M_in1,T_in1,'r');

%%
[gridM2,gridAoA2] =  ndgrid(plotM,plotT);



% Run First Stage =========================================================
const_firststage = 1;
addpath('../../FirstStage')
% addpath('../../DIDO_7.3.7')
% run startup.m
[FirstStageStates] = FirstStageProblem(alt(1),gamma(1),lat(1),zeta(1),const_firststage);
cd('../SecondStage/Combined - 2ns Stage Ascent and Return')
dlmwrite('FirstStage.txt', FirstStageStates);
copyfile('FirstStage.txt',sprintf('../ArchivedResults/%s/firststage_%s.txt',Timestamp,Timestamp))


%% Latitude Plot
figure(250)
plot(FirstStageStates(:,9))
plot(phi)
plot(ThirdStagePhi)
title('Latitude')

%% SAVE FIGS
saveas(figure(301),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'ThirdStage.fig']);
saveas(figure(101),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'FirstStage.fig']);
%%

% =========================================================================
% Troubleshooting Procedure
% =========================================================================

% 1: Check that you have posed your problem correctly ie. it is physically
% feasible and the bounds allow for a solution
% 2: Check for NaN values (check derivatives in Dynamics file while running)
% 3: Check guess, is it reasonable? Is it too close to the expected
% solution? Both can cause errors! Sometimes there is no real rhyme or
% reason to picking the correct guess, but a close bound to
% the expected solution has worked the most in my experience
% 4: Play with the no. of nodes, try both even and odd values
% 5: Play with scaling
% 6: Try all of the above in various combinations until it works!





