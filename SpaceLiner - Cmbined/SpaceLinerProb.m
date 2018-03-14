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

%% Get population interpolation
load PopInterp

auxdata.PopInterp = PopInterp;

%% Import Vehicle Config Data %%============================

Stage1.A = 93.25; %Reference Area in m² 
Stage1.mStruct = 205245.8;

% mixture ratio 6
% Stage1.T_SL = 1961*9*1e3; %for mixture ratio 6 with all boosters active
% Stage1.Isp_SL = 389;
% 
% Stage1.T_vac = 2206*9*1e3;
% Stage1.Isp_vac = 437;

% mixture ration 6.5
Stage1.T_SL = 2111*9*1e3; %for mixture ratio 6 with all boosters active
Stage1.Isp_SL = 390;

Stage1.T_vac = 2356*9*1e3;
Stage1.Isp_vac = 435;

Stage2.A = 461; %Reference Area in m² 
Stage2.mStruct = 134361.1;

% mixture ratio 6
% Stage2.T_SL = 1830*2*1e3; %for mixture ratio 6 with all boosters active
% Stage2.Isp_SL = 363;
% 
% Stage2.T_vac = 2268*2*1e3;
% Stage2.Isp_vac = 449;

% mixture ration 6.5
Stage2.T_SL = 1986*2*1e3; %for mixture ratio 6 with all boosters active
Stage2.Isp_SL = 367;

Stage2.T_vac = 2425*2*1e3;
Stage2.Isp_vac = 448;

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



Stage2_Cl_scattered = scatteredInterpolant(Stage2_MList_Full,Stage2_aoaList_Full,Stage2_Cl);
Stage2_Cd_scattered = scatteredInterpolant(Stage2_MList_Full,Stage2_aoaList_Full,Stage2_Cd);


Stage2_MList_interp = [Stage2_MList(1):0.1:Stage2_MList(end)];

[Stage2_MGrid,Stage2_aoaGrid] = meshgrid(Stage2_MList_interp,Stage2_aoaList);

for i = 1:length(Stage2_MList_interp)
    for j = 1:length(Stage2_aoaList)
        Stage2_Cl_Grid(j,i) = Stage2_Cl_scattered(Stage2_MList_interp(i),Stage2_aoaList(j));
        Stage2_Cd_Grid(j,i) = Stage2_Cd_scattered(Stage2_MList_interp(i),Stage2_aoaList(j));
    end
end

%% Create Gridded Interpolants

auxdata.interp.Stage1.Cl_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cl_Grid','linear');
auxdata.interp.Stage1.Cd_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cd_Grid','linear');

auxdata.interp.Stage2.Cl_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cl_Grid','linear');
auxdata.interp.Stage2.Cd_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cd_Grid','linear');


% auxdata.interp.Stage1.Cl_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cl_Grid','spline');
% auxdata.interp.Stage1.Cd_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cd_Grid','spline');
% 
% auxdata.interp.Stage2.Cl_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cl_Grid','spline');
% auxdata.interp.Stage2.Cd_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cd_Grid','spline');


%% Import Bounds %%========================================================
lonMin = -2*pi;         lonMax = -lonMin;

% latMin = -70*pi/180;  latMax = -latMin; %for japan launch

% latMin = -pi/2+0.0000001;  
% latMax = pi/2-0.0000001;

latMin = deg2rad(-89.5);  
latMax = deg2rad(89.5);

% lat0 = deg2rad(45);
% lon0 = deg2rad(145);


mission = 3


% lat0 = deg2rad(-23.3791); % Rockhampton
% lon0 = deg2rad(150.5100); % Rockhampton

if mission == 2
% lat0 = deg2rad(53.77); % Germany
% lon0 = deg2rad(8.6359);% Germany
% lat0 = deg2rad(54.661308); % Germany over water
% lon0 = deg2rad(6.847264);% Germany
lat0 = deg2rad(54.330731); % Germany over water
lon0 = deg2rad(7.763275);% Germany

%     lat0 = deg2rad(29.56043); %just west of florida (north)
% lon0 = deg2rad(-85.614341);
end

if mission == 1
lat0 = deg2rad(38.745095); %north korea
lon0 = deg2rad(128.272375);

% lat0 = deg2rad(36.065504); %south korea
% lon0 = deg2rad(129.665407);
end

if mission == 3
%     lat0 = deg2rad(27.219594); %just west of florida
% lon0 = deg2rad(-82.742209);

    lat0 = deg2rad(29.56043); %just west of florida (north)
lon0 = deg2rad(-85.614341);

%     lat0 = deg2rad(24.429926); %just west of florida (south)
% lon0 = deg2rad(-82.292144);

% lat0 = deg2rad(38.428333); %lisbon test
% lon0 = deg2rad(-9.706724);

% lat0 = deg2rad(54.330731); % Germany over water test
% lon0 = deg2rad(7.763275);% Germany

% lat0 = deg2rad(4.407410); %  over water test
% lon0 = deg2rad(-81.249230);% 

% lat0 = deg2rad(-38.95826); %  over water south test
% lon0 = deg2rad(-75.636537);% 

end

if mission == 4
   lat0 = deg2rad(-23.3791); % Rockhampton
lon0 = deg2rad(150.5100); % Rockhampton 
end

if mission == 5
   lat0 = deg2rad(-23.3791); % Rockhampton
lon0 = deg2rad(150.5100); % Rockhampton 
end

if mission == 6
   lat0 = deg2rad(-37.341524); % near mt gambier
lon0 = deg2rad(139.766669); % 
end

% lat0 = deg2rad(-23.3791); % Rockhampton
% lon0 = deg2rad(150.5100); % Rockhampton

%  lat0 = deg2rad(37.533151); % Close to Suzu, middle north Japan
% lon0 = deg2rad(137.276039);


% latF = deg2rad(53.77); % Germany
% lonF = deg2rad(8.6359);% Germany

if mission ==1 || mission ==4
latF = deg2rad(54.610440); % Germany over water
lonF = deg2rad(7.440526);% Germany

end

if mission == 2
%  latF = deg2rad(37.533151); % Close to Suzu, middle north Japan
% lonF = deg2rad(137.276039);
%  latF = deg2rad(29.823758); % china
% lonF = deg2rad(124.779639);

latF = deg2rad(38.745095); %north korea
lonF = deg2rad(128.272375);

%    latF = deg2rad(-23.3791); % Rockhampton test
% lonF = deg2rad(150.5100); % Rockhampton 

%    latF = deg2rad(16.697713); % manilla test
% lonF = deg2rad(123.120771); % 


%    latF = deg2rad(43.439108); % eastern tip of sapporo
% lonF = deg2rad(146.187017); % 


end

if mission == 3
   latF = deg2rad(-23.3791); % Rockhampton
lonF = deg2rad(150.5100); % Rockhampton 
end

if mission == 5
%     latF = deg2rad(27.219594); %just west of florida
% lonF = deg2rad(-82.742209);
    latF = deg2rad(29.56043); %just west of florida
lonF = deg2rad(-85.614341);
 
end

if mission ==6 
%     latF = deg2rad(38.428333); %lisbon
% lonF = deg2rad(-9.706724);

    latF = deg2rad(-4.973436); %atlantic
lonF = deg2rad(-18.291222);

end

auxdata.lon0 = lon0;

aoaMin = 0*pi/180;  aoaMax = 30*pi/180;

% bankMin_ascent = -10*pi/180; bankMax_ascent =   10*pi/180;
bankMin_ascent = -10*pi/180; bankMax_ascent =  10*pi/180;
bankMin_descent = -50*pi/180; bankMax_descent =   50*pi/180;

% Primal Bounds
if mission == 4 || mission == 5 || mission == 6
bounds.phase(1).state.lower = [0, lonMin, latMin, 0, -deg2rad(89), -pi, 0, aoaMin, bankMin_ascent];
bounds.phase(1).state.upper = [200000, lonMax, latMax, 10000, deg2rad(89), pi, 1.6038e+06, aoaMax, bankMax_ascent];
end
if mission ==1 
bounds.phase(1).state.lower = [0, lonMin, 0, 0, -deg2rad(89), -pi, 0, aoaMin, bankMin_ascent];
bounds.phase(1).state.upper = [200000, lonMax, latMax, 10000, deg2rad(89), pi, 1.6038e+06, aoaMax, bankMax_ascent];
end
if mission == 2
bounds.phase(1).state.lower = [0, -0.5, lat0-0.5, 0, -deg2rad(89), deg2rad(90), 0, aoaMin, bankMin_ascent];
bounds.phase(1).state.upper = [200000, 0.5, lat0+0.5, 10000, deg2rad(89), deg2rad(270), 1.6038e+06, aoaMax, bankMax_ascent];
end
if mission ==3 
bounds.phase(1).state.lower = [0, -1.5, lat0-1.5, 0, -deg2rad(89), 0, 0, aoaMin, bankMin_ascent];
bounds.phase(1).state.upper = [200000, 1.5, lat0+1.5, 10000, deg2rad(89), 2*pi, 1.6038e+06, aoaMax, bankMax_ascent];
end
% Initial States
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 95, deg2rad(80), deg2rad(70), 1.5038e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 105, deg2rad(89), deg2rad(80), 1.6038e+06, aoaMax];

% %Japan-Germany
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 100, deg2rad(80), deg2rad(60), 1.28e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 120, deg2rad(85), deg2rad(85), 1.32e+06, aoaMax];

%NKorea-Germany
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 100, deg2rad(80), deg2rad(60), 1.28e+06, aoaMin, bankMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 120, deg2rad(85), deg2rad(85), 1.32e+06, aoaMax, bankMax];

if mission ==1 || mission == 5|| mission == 6

bounds.phase(1).initialstate.lower = [500,-0.005, lat0-0.005, 50, deg2rad(87), deg2rad(-89), 213213, aoaMin, bankMin_ascent] ;
bounds.phase(1).initialstate.upper = [550,0.005, lat0+0.005, 60, deg2rad(88), deg2rad(89), 1.37e+06, aoaMax, bankMax_ascent];


% bounds.phase(1).initialstate.lower = [5900,-0.005, lat0-0.005, 250, deg2rad(63), deg2rad(-89), 213213, aoaMin, bankMin_ascent] ;
% bounds.phase(1).initialstate.upper = [6100,0.005, lat0+0.005, 260, deg2rad(64), deg2rad(89), 1.14e+06, aoaMax, bankMax_ascent];

end

if mission == 4 
bounds.phase(1).initialstate.lower = [500,0, lat0, 50, deg2rad(87), deg2rad(75), 213213, aoaMin, bankMin_ascent] ;
bounds.phase(1).initialstate.upper = [550,0, lat0, 60, deg2rad(88), deg2rad(89), 1.37e+06, aoaMax, bankMax_ascent];

end


if mission == 2
bounds.phase(1).initialstate.lower = [500,-0.005, lat0-0.005, 50, deg2rad(87), deg2rad(91), 213213, aoaMin, bankMin_ascent] ;
bounds.phase(1).initialstate.upper = [550,0.005, lat0+0.005, 60, deg2rad(88), deg2rad(269), 1.37e+06, aoaMax, bankMax_ascent];

end

if mission ==3 
bounds.phase(1).initialstate.lower = [500,-0.005, lat0-0.005, 50, deg2rad(86), 0, 1.30e+06, aoaMin, bankMin_ascent] ;
bounds.phase(1).initialstate.upper = [550,0.005, lat0+0.005, 70, deg2rad(88), 2*pi, 1.37e+06, aoaMax, bankMax_ascent];


% bounds.phase(1).initialstate.lower = [5900,-0.005, lat0-0.005, 250, deg2rad(63), deg2rad(91), 213213, aoaMin, bankMin_ascent] ;
% bounds.phase(1).initialstate.upper = [6100,0.005, lat0+0.005, 260, deg2rad(64), deg2rad(269), 1.14e+06, aoaMax, bankMax_ascent];

end

% bounds.phase(1).initialstate.lower = [1000,0, lat0, 100, deg2rad(80), deg2rad(-80), 1.28e+06, aoaMin, bankMin_ascent] ;
% bounds.phase(1).initialstate.upper = [1200,0, lat0, 120, deg2rad(85), deg2rad(85), 1.32e+06, aoaMax, bankMax_ascent];


%Germany-Japan
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 100, deg2rad(80), deg2rad(110), 1.28e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 120, deg2rad(85), deg2rad(130), 1.32e+06, aoaMax];


bounds.phase(1).finalstate.lower = bounds.phase(1).state.lower;
bounds.phase(1).finalstate.upper = bounds.phase(1).state.upper;


bounds.phase(2).initialstate.lower = bounds.phase(1).state.lower;
bounds.phase(2).initialstate.upper = bounds.phase(1).state.upper;

% Control Bounds
bounds.phase(1).control.lower = [deg2rad(-1), deg2rad(-2)];
bounds.phase(1).control.upper = [deg2rad(1), deg2rad(2)];

% Time Bounds
bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower = 0;
bounds.phase(1).finaltime.upper = 1000;

%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.

% bounds.phase(1).path.lower = [0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
% bounds.phase(1).path.upper = [60000, 2.5*9.81];



if mission == 1 || mission == 4 ||  mission == 6
% bounds.phase(1).path.lower = [0, 0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
% bounds.phase(1).path.upper = [60000, 7e6, 2.5*9.81];
bounds.phase(1).path.lower = [0, 0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
bounds.phase(1).path.upper = [40000, 1.3e6, 2.5*9.81];
end

if mission == 3 
bounds.phase(1).path.lower = [0, 0, 0]; % if using total acceleration, this might nee dto be in gs
bounds.phase(1).path.upper = [40000, 1.3e6, 2.5*9.81];

% bounds.phase(1).path.lower = [0, 0]; % if using total acceleration, this might nee dto be in gs
% bounds.phase(1).path.upper = [40000, 2.5*9.81];

end

if mission == 5 
bounds.phase(1).path.lower = [0, 0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
bounds.phase(1).path.upper = [40000, 1.3e6, 2.5*9.81];
end

if mission == 2
bounds.phase(1).path.lower = [0, 0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
bounds.phase(1).path.upper = [40000, 1.3e6, 2.5*9.81];
end

% bounds.phase(1).path.lower = [0, 0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
% bounds.phase(1).path.upper = [60000, 1.3e6, 2.5*9.81];

%% Bound integral if necessary
bounds.phase(1).integral.lower = 0;
bounds.phase(1).integral.upper = 1e9;


%% Set all phase bounds
bounds.phase(2).state = bounds.phase(1).state;

% if mission == 2
% bounds.phase(2).state.upper(3) = latMax;
% end

bounds.phase(2).finalstate = bounds.phase(1).finalstate;

bounds.phase(2).initialtime.lower = 0;
bounds.phase(2).initialtime.upper = 1000;

bounds.phase(2).finaltime = bounds.phase(1).finaltime;

bounds.phase(2).control = bounds.phase(1).control;

bounds.phase(2).path = bounds.phase(1).path;
% bounds.phase(2).path.lower(4) = -20;
% bounds.phase(2).path.upper(4) = 0;

bounds.phase(2).integral = bounds.phase(1).integral;

bounds.phase(3) = bounds.phase(2);
bounds.phase(4) = bounds.phase(2);
bounds.phase(5) = bounds.phase(2);
bounds.phase(6) = bounds.phase(2);
bounds.phase(7) = bounds.phase(2);
% bounds.phase(7) = bounds.phase(2);

% bounds.phase(7).finalstate.lower(7) = 213213; %set max propellant of orbiter at separation 
% bounds.phase(7).finalstate.upper(7) = 213213; %set max propellant of orbiter at separation 



bounds.phase(8).initialtime = bounds.phase(2).initialtime;
bounds.phase(8).finaltime = bounds.phase(2).finaltime;
bounds.phase(8).state = bounds.phase(2).state;
bounds.phase(8).initialstate = bounds.phase(2).initialstate;
bounds.phase(8).path = bounds.phase(2).path;
bounds.phase(8).integral = bounds.phase(2).integral;

% bounds.phase(8).initialstate.lower(7) = 213212; %set max propellant of orbiter at separation 
bounds.phase(8).initialstate.upper(7) = 213213; %set max propellant of orbiter at separation 



bounds.phase(8).finalstate = bounds.phase(2).finalstate;

% bounds.phase(8).finalstate.lower(1) = 100000;
% bounds.phase(8).finalstate.upper(1) = 100000;

bounds.phase(8).control = bounds.phase(2).control;


bounds.phase(9) = bounds.phase(8);


% bounds.phase(9).path.lower = [0, 0, 0]; % if using total acceleration, this might nee dto be in gs
% bounds.phase(9).path.upper = [40000, 1.3e6, 2.5*9.81];

% bounds.phase(9).path = bounds.phase(1).path;

bounds.phase(9).state.lower(9) = bankMin_descent;
bounds.phase(9).state.upper(9) = bankMax_descent;
bounds.phase(9).initialtime.upper = 1000;
bounds.phase(9).finaltime.upper = 10000;


bounds.phase(9).state.lower(2) = lonMin;
bounds.phase(9).state.lower(3) = latMin;
bounds.phase(9).state.upper(2) = lonMax;
bounds.phase(9).state.upper(3) = latMax;


% if mission == 3 
% bounds.phase(9).state.lower(2) = lonMin;
% bounds.phase(9).state.lower(3) = -0.5;
% bounds.phase(9).state.upper(2) = lonMax;
% bounds.phase(9).state.upper(3) = 0.5;
% end

% End States

% bounds.phase(8).finalstate.lower = [70000, lonMin, latMin, 0, 0, -2*pi, 0, aoaMin];
% bounds.phase(8).finalstate.upper = [150000, lonMax, latMax, 15000, deg2rad(80), 2*pi,1.6038e+06, aoaMax];
% 
% bounds.phase(8).finalstate.lower = [70000, lonMin, latMin, 5000, 0, -2*pi, 0, aoaMin];
% bounds.phase(8).finalstate.upper = [70000, lonMax, latMax, 15000, deg2rad(0), 2*pi,1.6038e+06, aoaMax];

% bounds.phase(8).finalstate.lower = [70000, lonMin, latMin, 7000, 0,
% -2*pi, 0, aoaMin]; % worked 
% bounds.phase(8).finalstate.upper = [70000, lonMax, latMax, 7000, deg2rad(0), 2*pi,1.6038e+06, aoaMax];

% Japan-Germany 
% bounds.phase(8).finalstate.lower = [69000, lonMin, latMin, 6800,deg2rad(-1), deg2rad(75), 0, aoaMin]; % Japan-Germany
% bounds.phase(8).finalstate.upper = [71000, lonMax, latMax, 7200, deg2rad(1), deg2rad(75),1.6038e+06, aoaMax];

% Germany- Japan
% bounds.phase(8).finalstate.lower = [69000, lonMin, latMin, 6800, deg2rad(-1), deg2rad(120), 0, aoaMin];
% bounds.phase(8).finalstate.upper = [71000, lonMax, latMax, 7200, deg2rad(1), deg2rad(120),1.6038e+06, aoaMax];

% korea-Germany 
% bounds.phase(9).finalstate.lower = [99000, lonMin, latMin, 6800,deg2rad(-30), deg2rad(75), 0, aoaMin, bankMin]; % Japan-Germany
% bounds.phase(9).finalstate.upper = [101000, lonMax, latMax, 7200, deg2rad(30), deg2rad(75),1.6038e+06, aoaMax, bankMax];

% bounds.phase(9).finalstate.lower = [1000, lonMin, latMin, 100,deg2rad(-89), deg2rad(10), 0, aoaMin, bankMin_descent]; % Japan-Germany
% bounds.phase(9).finalstate.upper = [100000, lonMax, latMax, 10000, deg2rad(89), deg2rad(10),1.6038e+06, aoaMax, bankMax_descent];

if mission == 1|| mission == 4 || mission == 6

bounds.phase(9).finalstate.lower = [1000, lonF-lon0+2*pi-0.005, latF-0.005, 100,deg2rad(-20), -deg2rad(90), 0, aoaMin, bankMin_descent]; % Japan-Germany
bounds.phase(9).finalstate.upper = [5000, lonF-lon0+2*pi+0.005, latF+0.005, 200, deg2rad(20), deg2rad(90),213213, aoaMax, bankMax_descent];

end

% if mission == 2
% 
% bounds.phase(9).finalstate.lower = [1000, lonF-lon0-2*pi, latF, 100,deg2rad(-20), -deg2rad(90), 0, aoaMin, bankMin_descent]; % Japan-Germany
% bounds.phase(9).finalstate.upper = [5000, lonF-lon0-2*pi, latF, 200, deg2rad(20), deg2rad(90),1.6038e+06, aoaMax, bankMax_descent];
% end

if mission == 2

bounds.phase(9).finalstate.lower = [1000, lonF-lon0-2*pi, latF, 100,deg2rad(-20), deg2rad(90), 0, aoaMin, bankMin_descent]; % Japan-Germany
bounds.phase(9).finalstate.upper = [5000, lonF-lon0-2*pi, latF, 200, deg2rad(20), deg2rad(270),213213, aoaMax, bankMax_descent];
end

if mission == 3 

bounds.phase(9).finalstate.lower = [1000, lonF-lon0-2*pi-0.005, latF-0.005, 100,deg2rad(-20), deg2rad(90), 0, aoaMin, bankMin_descent]; % FLorida-Aus
bounds.phase(9).finalstate.upper = [5000, lonF-lon0-2*pi+0.005, latF+0.005, 200, deg2rad(20), deg2rad(270),213213, aoaMax, bankMax_descent];
end

if mission == 5

bounds.phase(9).finalstate.lower = [1000, lonF-lon0+2*pi-0.005, latF-0.005, 100,deg2rad(-20), deg2rad(-90), 0, aoaMin, bankMin_descent]; 
bounds.phase(9).finalstate.upper = [5000, lonF-lon0+2*pi+0.005, latF+0.005, 200, deg2rad(20), deg2rad(90),213213, aoaMax, bankMax_descent];
end

% bounds.phase(9).finalstate.lower = [99000, lonMin, latMin, 6800,deg2rad(-1), deg2rad(60), 0, aoaMin, bankMin_descent]; % Japan-Germany
% bounds.phase(9).finalstate.upper = [101000, lonMax, latMax, 7200, deg2rad(1), deg2rad(60),1.6038e+06, aoaMax, bankMax_descent];


% bounds.phase(9).finalstate.lower = bounds.phase(1).state.lower;
% bounds.phase(9).finalstate.upper = bounds.phase(1).state.upper;

% bounds.phase(9).finalstate.lower(2) =lonF-lon0+2*pi;
% bounds.phase(9).finalstate.upper(2) =lonF-lon0+2*pi;
% 
% bounds.phase(9).finalstate.lower(3) =latF;
% bounds.phase(9).finalstate.upper(3) =latF;
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

bounds.eventgroup(7).lower = zeros(1,9);
bounds.eventgroup(7).upper = zeros(1,9);

bounds.eventgroup(8).lower = zeros(1,10);
bounds.eventgroup(8).upper = zeros(1,10);

bounds.eventgroup(9).lower = [ones(1,8) 1];
bounds.eventgroup(9).upper = [1000*ones(1,8) 10000];

%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(1).state(:,1)   = [400;10000];

% guess.phase(1).state(:,2)   = [2.5;2.55];
% guess.phase(1).state(:,3)   = [0.7;.75]; %japan - germany

% guess.phase(1).state(:,2)   = [2.238;2.3];
% guess.phase(1).state(:,2)   = [0;0.1];

if mission == 1|| mission == 4|| mission == 5  
    
guess.phase(1).state(:,2)   = [0;0.01];
end

if mission == 2 || mission == 3 || mission == 6  
    
guess.phase(1).state(:,2)   = [0;-0.01];
end

% guess.phase(1).state(:,3)   = [0.67;.75]; %japan - germany

if mission == 1 || mission == 5 || mission == 6  
guess.phase(1).state(:,3)   = [lat0;lat0+0.1]; % rockhampton-germany
end

if mission == 2 || mission == 3
guess.phase(1).state(:,3)   = [lat0;lat0+0.0]; 
end

% guess.phase(1).state(:,2)   = [0.15;0.1]; %germany-japan
% guess.phase(1).state(:,3)   = [0.95;.1];


guess.phase(1).state(:,4)   = [100,1100];
guess.phase(1).state(:,5)   = [deg2rad(85),deg2rad(85)];

if mission == 1 || mission == 4 
guess.phase(1).state(:,6)   = [deg2rad(75),deg2rad(75)]; %Japan-germany
% guess.phase(1).state(:,6)   = [deg2rad(15),deg2rad(15)]; %Japan-germany

end
if mission == 2
guess.phase(1).state(:,6)   = [deg2rad(110),deg2rad(110)]; %Germany-Japan
end

if mission == 3
guess.phase(1).state(:,6)   = [deg2rad(180),deg2rad(180)];
end

if mission == 5 
guess.phase(1).state(:,6)   = [deg2rad(-10),deg2rad(-10)];
end
if  mission == 6  
guess.phase(1).state(:,6)   = [deg2rad(-80),deg2rad(-80)];
end


guess.phase(1).state(:,7) 	= [1.37e+06, 1.17e+06];
guess.phase(1).state(:,8)   = [1*pi/180; 10*pi/180];
guess.phase(1).state(:,9)   = [deg2rad(0);deg2rad(0)];


guess.phase(1).control      = [[0;0],[0;0]];
if mission ==3
guess.phase(1).time          = [0;50];
else
guess.phase(1).time          = [0;150];
end

guess.phase(1).integral = 0;
% Tire stages together
% bounds.eventgroup(1).lower = [zeros(1,10)];
% bounds.eventgroup(1).upper = [zeros(1,10)]; 

if mission == 1
guess.phase(2).state = guess.phase(1).state + [15000 15000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(3).state = guess.phase(2).state+ [15000 15000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(4).state = guess.phase(3).state+ [15000 15000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(5).state = guess.phase(4).state+ [15000 15000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(6).state = guess.phase(5).state+ [15000 15000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(7).state = guess.phase(6).state+ [15000 15000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(8).state = guess.phase(7).state+ [15000 15000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(9).state = guess.phase(8).state+ [15000 -100000; 0.05 lonF-lon0+2*pi;.35 .05; 1000 -5000; -deg2rad(0) -deg2rad(60);0 -deg2rad(150);0 0;deg2rad(10) deg2rad(10);deg2rad(30) deg2rad(30)]';
end

if mission == 2
guess.phase(2).state = guess.phase(1).state + [8000 10000; -0.0 -0.0;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-200000 -200000;0 0;0 0]';
guess.phase(3).state = guess.phase(2).state+ [10000 10000; -0.0 -0.0;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-200000 -200000;0 0;0 0]';
guess.phase(4).state = guess.phase(3).state+ [10000 10000;-0.0 -0.0;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-200000 -200000;0 0;0 0]';
guess.phase(5).state = guess.phase(4).state+ [10000 10000; -0.0 -0.;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-200000 -200000;0 0;0 0]';
guess.phase(6).state = guess.phase(5).state+ [10000 10000; -0.0 -0.0;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-200000 -200000;0 0;0 0]';
guess.phase(7).state = guess.phase(6).state+ [10000 10000; -0.0 -0.0;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-200000 -200000;0 0;0 0]';
guess.phase(8).state = guess.phase(7).state+ [10000 10000; -0.0 -0.0;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-200000 -200000;0 0;0 0]';
guess.phase(9).state = guess.phase(8).state+ [10000 -70000; -0.0 lonF-lon0-2*pi;0 -.7; 1000 -5000; -deg2rad(5) -deg2rad(30);0 deg2rad(70);0 0;deg2rad(10) deg2rad(10);deg2rad(-30) deg2rad(-30)]';
end

if mission == 3
guess.phase(2).state = guess.phase(1).state + [8000 10000; -0.01 -0.01;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-100000 -100000;0 0;0 0]';
guess.phase(3).state = guess.phase(2).state+ [10000 10000; -0.01 -0.01;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-100000 -100000;0 0;0 0]';
guess.phase(4).state = guess.phase(3).state+ [10000 10000;-0.01 -0.01;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-100000 -100000;0 0;0 0]';
guess.phase(5).state = guess.phase(4).state+ [10000 10000; -0.01 -0.01;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-100000 -100000;0 0;0 0]';
guess.phase(6).state = guess.phase(5).state+ [10000 10000; -0.01 -0.01;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-100000 -100000;0 0;0 0]';
guess.phase(7).state = guess.phase(6).state+ [10000 10000; -0.01 -0.01;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-100000 -100000;0 0;0 0]';
guess.phase(8).state = guess.phase(7).state+ [10000 10000; -0.01 -0.01;.0 .0; 1000 1000; -deg2rad(5) -deg2rad(5);0 0;-100000 -100000;0 0;0 0]';
guess.phase(9).state = guess.phase(8).state+ [10000 -60000; -0.01 lonF-lon0-2*pi;lat0 latF; 1000 -5000; -deg2rad(5) -deg2rad(30);0 deg2rad(0);0 0;deg2rad(10) deg2rad(10);0 0]';
end

if mission == 4 
guess.phase(2).state = guess.phase(1).state + [10000 10000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(3).state = guess.phase(2).state+ [10000 10000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(4).state = guess.phase(3).state+ [10000 10000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(5).state = guess.phase(4).state+ [10000 10000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(6).state = guess.phase(5).state+ [10000 10000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(7).state = guess.phase(6).state+ [10000 10000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(8).state = guess.phase(7).state+ [10000 10000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(9).state = guess.phase(8).state+ [10000 -60000; 0.05 lonF-lon0+2*pi;1.5 latF; 1000 -5000; -deg2rad(0) -deg2rad(60);0 -deg2rad(150);0 0;deg2rad(10) deg2rad(10);deg2rad(30) deg2rad(30)]';
end

if mission == 5 
guess.phase(2).state = guess.phase(1).state + [10000 10000; 0.05 0.05;.0 .0; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(3).state = guess.phase(2).state+ [10000 10000; 0.05 0.05;.0 .0; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(4).state = guess.phase(3).state+ [10000 10000; 0.05 0.05;.0 .0; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(5).state = guess.phase(4).state+ [10000 10000; 0.05 0.05;.0 .0; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(6).state = guess.phase(5).state+ [10000 10000; 0.05 0.05;.0 .0; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(7).state = guess.phase(6).state+ [10000 10000; 0.05 0.05;.0 .0; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(8).state = guess.phase(7).state+ [10000 10000; 0.05 0.05;.0 .0; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(9).state = guess.phase(8).state+ [10000 -60000; 0.05 lonF-lon0+2*pi;0 latF; 1000 -5000; -deg2rad(0) -deg2rad(60);0 deg2rad(0);0 0;deg2rad(10) deg2rad(10);deg2rad(0) deg2rad(0)]';
end

if mission == 6 
guess.phase(2).state = guess.phase(1).state + [15000 15000; -0.05 -0.05;-.05 -.05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(3).state = guess.phase(2).state+ [15000 15000; -0.05 -0.05;-.05 -.05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(4).state = guess.phase(3).state+ [15000 15000;-0.05 -0.05;-.05 -.05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(5).state = guess.phase(4).state+ [15000 15000; -0.05 -0.05;-.05 -.05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(6).state = guess.phase(5).state+ [15000 15000; -0.05 -0.05;-.05 -.05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(7).state = guess.phase(6).state+ [15000 15000; -0.05 -0.05;-.05 -.05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(8).state = guess.phase(7).state+ [15000 15000; -0.05 -0.05;-.05 -.05; 1000 1000; -deg2rad(0) -deg2rad(0);0 0;-200000 -200000;0 0;0 0]';
guess.phase(9).state = guess.phase(8).state+ [15000 -70000; -0.05 lonF-lon0+2*pi;-1.35 -1.05; 1000 -5000; -deg2rad(0) -deg2rad(60);0 deg2rad(70);0 0;deg2rad(10) deg2rad(10);deg2rad(-30) deg2rad(-30)]';
end

guess.phase(2).control = guess.phase(1).control;
guess.phase(3).control = guess.phase(2).control;
guess.phase(4).control = guess.phase(3).control;
guess.phase(5).control = guess.phase(4).control;
guess.phase(6).control = guess.phase(5).control;
guess.phase(7).control = guess.phase(6).control;
guess.phase(8).control = guess.phase(7).control;
guess.phase(9).control = guess.phase(8).control;

if mission ==3
guess.phase(2).time = guess.phase(1).time + 10; 
guess.phase(3).time = guess.phase(2).time + 10;
guess.phase(4).time = guess.phase(3).time + 10;
guess.phase(5).time = guess.phase(4).time + 10;
guess.phase(6).time = guess.phase(5).time + 10;
guess.phase(7).time = guess.phase(6).time + 10;
guess.phase(8).time = guess.phase(7).time + 200;
guess.phase(9).time = guess.phase(8).time + 3000;
else
    guess.phase(2).time = guess.phase(1).time + 30; 
guess.phase(3).time = guess.phase(2).time + 30;
guess.phase(4).time = guess.phase(3).time + 30;
guess.phase(5).time = guess.phase(4).time + 30;
guess.phase(6).time = guess.phase(5).time + 30;
guess.phase(7).time = guess.phase(6).time + 30;
guess.phase(8).time = guess.phase(7).time + 30;
guess.phase(9).time = guess.phase(8).time + 3000;
end

guess.phase(2).integral = guess.phase(1).integral;
guess.phase(3).integral = guess.phase(2).integral;
guess.phase(4).integral = guess.phase(3).integral;
guess.phase(5).integral = guess.phase(4).integral;
guess.phase(6).integral = guess.phase(5).integral;
guess.phase(7).integral = guess.phase(6).integral;
guess.phase(8).integral = guess.phase(7).integral;
guess.phase(9).integral = guess.phase(8).integral;

%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 3;
mesh.colpointsmin = 3;
mesh.colpointsmax = 250;
mesh.tolerance    = 1e-4;


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
setup.nlp.ipoptoptions.maxiterations = 700;
setup.derivatives.supplier           = 'sparseCD';
% setup.derivatives.derivativelevel    = 'second';
% setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.derivativelevel    = 'first';
% setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
% setup.scales.method                  = 'automatic-guessUpdate';
setup.scales.method                  = 'automatic-hybridUpdate';
setup.derivatives.dependencies      = 'full';
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

controls.Alpha = Alpha';




% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 

% throttle = ones(length(alt),1);
% 
% throttle(time>=output.result.solution.phase(1).time(end)) = 0.891;
% 
% throttle(time>=output.result.solution.phase(2).time(end)) = 0.812;
% 
% throttle(time>=output.result.solution.phase(3).time(end)) = .7333;
% 
% throttle(time>=output.result.solution.phase(4).time(end)) = .6545;
% 
% throttle(time>=output.result.solution.phase(5).time(end)) = .5757;
% 
% throttle(time>=output.result.solution.phase(6).time(end)) = 0.496;
% 
% throttle(time>=output.result.solution.phase(7).time(end)) = 1; %after separation
% 
% throttle(time>=output.result.solution.phase(8).time(end)) = 0; %after separation


% Watch that this is throttleing at the right time / being interpolated
% correctly in the forward sims
throttle = ones(length(alt),1);
throttle_temp = 1;

stage = ones(1,length(time));
% stage(time>=output.result.solution.phase(7).time(end)) = 2;
stage_temp = 1;
for i = 1:length(throttle)
    throttle(i) = throttle_temp;
    stage(i) = stage_temp;
    
    if time(i) == output.result.solution.phase(1).time(end)
        throttle_temp = 0.891;
        stage_temp = 1;
    end
    
    if time(i) == output.result.solution.phase(2).time(end)
        throttle_temp = 0.812;
        stage_temp = 1;
    end
    
    if time(i) == output.result.solution.phase(3).time(end)
        throttle_temp =.7333;
        stage_temp = 1;
    end
    
    if time(i) == output.result.solution.phase(4).time(end)
        throttle_temp = .6545;
        stage_temp = 1;
    end
    
    if time(i) == output.result.solution.phase(5).time(end)
        throttle_temp = .5757;
        stage_temp = 1;
    end
    
    if time(i) == output.result.solution.phase(6).time(end)
        throttle_temp = 0.496;
        stage_temp = 1;
    end
    
    if time(i) == output.result.solution.phase(7).time(end)
        throttle_temp = 1;
        stage_temp = 2;
    end
    
    if time(i) == output.result.solution.phase(8).time(end)
        throttle_temp = 0;
        stage_temp = 2;
    end
    

end





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

% stage = 1; %will be wrong for orbiter

stage_forward = ones(1,length(time_forward));
stage_forward(time_forward>=output.result.solution.phase(7).time(end)) = 2;


[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time_forward,Alpha_forward,f_t),ControlInterp(time_forward,eta_forward,f_t),ControlInterp(time_forward,throttle_forward,f_t),ControlInterp(time_forward,stage_forward,f_t)),0:time(end),forward0);
% [f_t, f_y] = ode23(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time,Alpha,f_t),ControlInterp(time,throttle,f_t),time(end)),time(1:end),forward0);


% dt = 1;
% f_y = forward0;
% temp = 1;
% for t_temp = 0:dt:time(end)
%     
%     df_y = VehicleModel_forward(t_temp, f_y(temp,:),auxdata,ControlInterp(time,Alpha,t_temp),ControlInterp(time,throttle,t_temp),time(end));
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
plot(time,alt);

subplot(7,1,3)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time,gamma);


subplot(7,1,4)
hold on
plot(f_t(1:end),f_y(:,3));
plot(time,v);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,4));
plot(time,zeta);

subplot(7,1,7)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time,mFuel);

figure()
hold on
plot(f_y(:,6),f_y(:,5));
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


T = [];
M = [];
Fd = [];
L = [];
q1 = [];
m = [];
gammadot = [];
azidot = [];
for i = 1:length(alt)
    phase_temp.state = phase.state(i,:);
[altdot,londot,latdot,gammadot(i),vdot(i),azidot(i), q(i), M(i), Fd(i), rho,L(i),Fueldt1,T(i),Isp1,Isp2,m(i),heating_rate(i),total_acceleration(i),Cl(i),Cd(i)] = SpaceLinerVehicleModel(time(i),phase_temp,throttle(i),auxdata,stage(i));

end


%
altpop     = alt;
lonpop     = lon;
latpop     = lat;

lonpop = lonpop + auxdata.lon0;
lonpop(lonpop > pi) = lonpop(lonpop > pi) - 2*pi;
lonpop(lonpop < -pi) = lonpop(lonpop < -pi) + 2*pi;

AltCost1 = (80000-altpop)/10000;
AltCost1(alt>80000) = 0;
% AltCost1(altpop>90000) = gaussmf(altpop(altpop>90000),[10000 90000]);

pop1 = auxdata.PopInterp(rad2deg(lonpop),rad2deg(latpop));

popCost1 = pop1.*AltCost1; % for flights which go over large amounts of




figure(201)
subplot(5,2,1)
hold on
plot(time,alt/1000)
xlabel('time (s)')
ylabel('altitude (km)')

subplot(5,2,2)
hold on
plot(v,alt/1000)
xlabel('velocity (m/s)')
ylabel('altitude (km)')


subplot(5,2,3)
hold on
plot(time,rad2deg(Alpha))
plot(time,rad2deg(gamma))
plot(time,rad2deg(eta))
legend('Angle of Attack','Trajectory Angle', 'Bank Angle')
xlabel('Time (s)')
ylabel('(deg)')

subplot(5,2,4)
hold on
plot(time,T/1000)
plot(time,Fd/1000)
plot(time,L/1000)
xlabel('time (s)')
ylabel('kN')
legend('Thrust','Drag','Lift')

subplot(5,2,5)
hold on
plot(time,total_acceleration/9.81)
xlabel('time (s)')
ylabel('acceleration (g)')

subplot(5,2,6)
hold on
plot(time,m)
xlabel('time (s)')
ylabel('Mass (kg)')

subplot(5,2,7)
hold on
plot(time,M)
xlabel('time (s)')
ylabel('Mach no.')

subplot(5,2,8)
hold on
plot(time,q)
xlabel('time (s)')
ylabel('Dynamic Pressure (kPa)')

subplot(5,2,9)
hold on
plot(time,heating_rate/1e6)
xlabel('time (s)')
ylabel('Heating Rate (MW/m^2)')

figure(230)
hold on
plot3(lon,lat,alt)
plot3(lon,lat,popCost1)


%%

    
    figure(231)
hold on

axesm('pcarree','Origin',[0 rad2deg(lon0) 0])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
% plotm(rad2deg(lat),rad2deg(lon+lon0))
plotm(rad2deg(lat),rad2deg(lon+lon0))
    
    cities = shaperead('worldcities', 'UseGeoCoords', true);
lats = extractfield(cities,'Lat');
lons = extractfield(cities,'Lon');
geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r',...
        'MarkerSize', 2)
% 
%% save file
contents = {'time', 'v', 'gamma', 'zeta', 'alt', 'lon', 'lat', 'mFuel', 'PlaceHolder', 'Thrust', 'PlaceHolder', 'Accel', 'q', 'heating_rate', 'PlaceHolder' , 'PlaceHolder', 'PlaceHolder', 'PlaceHolder', 'PlaceHolder', 'L', 'D', 'M','eta', 'Alpha' };
units = {'s', 'km/s', 'deg', 'deg', 'km', 'deg', 'deg', 'kg', '-', 'N', '-', 'm/s^2', 'pa', 'W/m^2', '-', '-', '-', '-', '-', 'N', 'N', 'Mach', 'deg', 'deg'};

result = [time' v'/1000 rad2deg(gamma)' rad2deg(zeta)' alt'/1000 rad2deg(lon'+lon0) rad2deg(lat') mFuel' 0.*mFuel' T' 0.*mFuel' total_acceleration' q' heating_rate' 0.*mFuel' 0.*mFuel' 0.*mFuel' 0.*mFuel' 0.*mFuel' L' Fd' Cl' Cd' M'  rad2deg(eta)' rad2deg(Alpha)' ];
    delete('out')
    
fid=fopen('out','wt');

for i=1:length(contents)
      fprintf(fid,'%s\t',contents{i})

end

for i=1:length(units)
    if i == 1
      fprintf(fid,'\n%s\t',units{i})
    else
        fprintf(fid,'%s\t',units{i})
    end

end

    
dlmwrite('out', result, '-append','delimiter', '\t','roffset',1);
    

    
