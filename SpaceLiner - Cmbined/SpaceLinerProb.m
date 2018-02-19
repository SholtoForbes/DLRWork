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

auxdata.interp.Stage1.Cl_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cl_Grid','linear');
auxdata.interp.Stage1.Cd_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cd_Grid','linear');

auxdata.interp.Stage2.Cl_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cl_Grid','linear');
auxdata.interp.Stage2.Cd_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cd_Grid','linear');


%% Import Bounds %%========================================================
lonMin = -pi;         lonMax = -lonMin;

% latMin = -70*pi/180;  latMax = -latMin; %for japan launch

latMin = -pi/2+0.0000001;  
latMax = pi/2-0.0000001;

% lat0 = deg2rad(45);
% lon0 = deg2rad(145);

% lat0 = deg2rad(-23.3791); % Rockhampton
% lon0 = deg2rad(150.5100); % Rockhampton

% lat0 = deg2rad(53.77); % Germany
% lon0 = deg2rad(8.6359);% Germany

lat0 = deg2rad(38.745095); %north korea
lon0 = deg2rad(128.272375);

%  lat0 = deg2rad(37.533151); % Close to Suzu, middle north Japan
% lon0 = deg2rad(137.276039);


latF = deg2rad(53.77); % Germany
lonF = deg2rad(8.6359);% Germany


aoaMin = 0;  aoaMax = 20*pi/180;
% bankMin1 = -50*pi/180; bankMax1 =   50*pi/180;

% Primal Bounds
bounds.phase(1).state.lower = [0, lonMin, latMin, 0, -deg2rad(50), -pi, 0, aoaMin];
bounds.phase(1).state.upper = [200000, lonMax, latMax, 15000, deg2rad(85), pi, 1.6038e+06, aoaMax];

% Initial States
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 95, deg2rad(80), deg2rad(70), 1.5038e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 105, deg2rad(89), deg2rad(80), 1.6038e+06, aoaMax];

% %Japan-Germany
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 100, deg2rad(80), deg2rad(60), 1.28e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 120, deg2rad(85), deg2rad(85), 1.32e+06, aoaMax];

%NKorea-Germany
bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 100, deg2rad(80), deg2rad(-10), 1.28e+06, aoaMin] ;
bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 120, deg2rad(85), deg2rad(10), 1.32e+06, aoaMax];


%Germany-Japan
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 100, deg2rad(80), deg2rad(110), 1.28e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 120, deg2rad(85), deg2rad(130), 1.32e+06, aoaMax];


bounds.phase(1).finalstate.lower = bounds.phase(1).state.lower;
bounds.phase(1).finalstate.upper = bounds.phase(1).state.upper;


bounds.phase(2).initialstate.lower = bounds.phase(1).state.lower;
bounds.phase(2).initialstate.upper = bounds.phase(1).state.upper;

% Control Bounds
bounds.phase(1).control.lower = [deg2rad(-1)];
bounds.phase(1).control.upper = [deg2rad(1)];

% Time Bounds
bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower = 50;
bounds.phase(1).finaltime.upper = 5000;

%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.

% bounds.phase(1).path.lower = [0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
% bounds.phase(1).path.upper = [40000, 2.5*9.81];

bounds.phase(1).path.lower = [0, -2.5*9.81]; % if using total acceleration, this might nee dto be in gs
bounds.phase(1).path.upper = [60000, 2.5*9.81];

%% Bound integral if necessary
bounds.phase(1).integral.lower = 0;
bounds.phase(1).integral.upper = 1e9;


%% Set all phase bounds
bounds.phase(2).state = bounds.phase(1).state;
bounds.phase(2).finalstate = bounds.phase(1).finalstate;

bounds.phase(2).initialtime.lower = 0;
bounds.phase(2).initialtime.upper = 5000;

bounds.phase(2).finaltime = bounds.phase(1).finaltime;

bounds.phase(2).control = bounds.phase(1).control;

bounds.phase(2).path = bounds.phase(1).path;

bounds.phase(2).integral = bounds.phase(1).integral;

bounds.phase(3) = bounds.phase(2);
bounds.phase(4) = bounds.phase(2);
bounds.phase(5) = bounds.phase(2);
bounds.phase(6) = bounds.phase(2);
bounds.phase(7) = bounds.phase(2);
bounds.phase(7) = bounds.phase(2);



bounds.phase(8).initialtime = bounds.phase(2).initialtime;
bounds.phase(8).finaltime = bounds.phase(2).finaltime;
bounds.phase(8).state = bounds.phase(2).state;
bounds.phase(8).initialstate = bounds.phase(2).initialstate;
bounds.phase(8).path = bounds.phase(2).path;
bounds.phase(8).integral = bounds.phase(2).integral;
bounds.phase(8).initialstate.upper(7) = 213213; %set max propellant of orbiter at separation 
bounds.phase(8).finalstate = bounds.phase(2).finalstate;

bounds.phase(8).control = bounds.phase(2).control;


bounds.phase(9) = bounds.phase(8);



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
bounds.phase(9).finalstate.lower = [1000, lon0+1.5, latMin, 100,deg2rad(-1), deg2rad(60), 0, aoaMin]; % Japan-Germany
bounds.phase(9).finalstate.upper = [101000, lon0+1.5, latMax, 1000, deg2rad(1), deg2rad(60),1.6038e+06, aoaMax];

%% Event bounds
bounds.eventgroup(1).lower = zeros(1,9);
bounds.eventgroup(1).upper = zeros(1,9);

bounds.eventgroup(2).lower = zeros(1,9);
bounds.eventgroup(2).upper = zeros(1,9);

bounds.eventgroup(3).lower = zeros(1,9);
bounds.eventgroup(3).upper = zeros(1,9);

bounds.eventgroup(4).lower = zeros(1,9);
bounds.eventgroup(4).upper = zeros(1,9);

bounds.eventgroup(5).lower = zeros(1,9);
bounds.eventgroup(5).upper = zeros(1,9);

bounds.eventgroup(6).lower = zeros(1,9);
bounds.eventgroup(6).upper = zeros(1,9);

bounds.eventgroup(7).lower = zeros(1,9);
bounds.eventgroup(7).upper = zeros(1,9);

bounds.eventgroup(8).lower = zeros(1,9);
bounds.eventgroup(8).upper = zeros(1,9);

bounds.eventgroup(9).lower = ones(1,9);
bounds.eventgroup(9).upper = 1000*ones(1,9);

%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(1).state(:,1)   = [2000;7000];

% guess.phase(1).state(:,2)   = [2.5;2.55];
% guess.phase(1).state(:,3)   = [0.7;.75]; %japan - germany

guess.phase(1).state(:,2)   = [2.238;2.3];
guess.phase(1).state(:,3)   = [0.67;.75]; %japan - germany

% guess.phase(1).state(:,2)   = [0.15;0.1]; %germany-japan
% guess.phase(1).state(:,3)   = [0.95;.1];


guess.phase(1).state(:,4)   = [100,1100];
guess.phase(1).state(:,5)   = [deg2rad(80),deg2rad(80)];

% guess.phase(1).state(:,6)   = [deg2rad(75),deg2rad(80)]; %Japan-germany
% guess.phase(1).state(:,6)   = [deg2rad(117),deg2rad(120)]; %Germany-Japan
% guess.phase(1).state(:,6)   = [deg2rad(71),deg2rad(80)]; %korea-germany
guess.phase(1).state(:,6)   = [deg2rad(1),deg2rad(1)]; %korea-germany

guess.phase(1).state(:,7) 	= [1.6038e+06, 1.500e+06];
guess.phase(1).state(:,8)   = [1*pi/180; 1*pi/180];
% guess.phase(1).state(:,9)   = [deg2rad(0);deg2rad(0)];


guess.phase(1).control      = [[0;0]];
guess.phase(1).time          = [0;50];


guess.phase(1).integral = 0;
% Tire stages together
% bounds.eventgroup(1).lower = [zeros(1,10)];
% bounds.eventgroup(1).upper = [zeros(1,10)]; 

guess.phase(2).state = guess.phase(1).state + [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';
guess.phase(3).state = guess.phase(2).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';
guess.phase(4).state = guess.phase(3).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';
guess.phase(5).state = guess.phase(4).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';
guess.phase(6).state = guess.phase(5).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';
guess.phase(7).state = guess.phase(6).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';
guess.phase(8).state = guess.phase(7).state+ [5000 5000; 0.05 0.05;.05 .05; 1000 1000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';
guess.phase(9).state = guess.phase(7).state+ [5000 -50000; 0.05 0.05;.05 .05; 1000 -2000; -deg2rad(10) -deg2rad(10);0 0;-100000 -100000;0 0]';


guess.phase(2).control = guess.phase(1).control;
guess.phase(3).control = guess.phase(2).control;
guess.phase(4).control = guess.phase(3).control;
guess.phase(5).control = guess.phase(4).control;
guess.phase(6).control = guess.phase(5).control;
guess.phase(7).control = guess.phase(6).control;
guess.phase(8).control = guess.phase(7).control;
guess.phase(9).control = guess.phase(7).control;

guess.phase(2).time = guess.phase(1).time + 5; 
guess.phase(3).time = guess.phase(2).time + 5;
guess.phase(4).time = guess.phase(3).time + 5;
guess.phase(5).time = guess.phase(4).time + 5;
guess.phase(6).time = guess.phase(5).time + 5;
guess.phase(7).time = guess.phase(6).time + 5;
guess.phase(8).time = guess.phase(7).time + 5;
guess.phase(9).time = guess.phase(7).time + 100;

guess.phase(2).integral = guess.phase(1).integral;
guess.phase(3).integral = guess.phase(2).integral;
guess.phase(4).integral = guess.phase(3).integral;
guess.phase(5).integral = guess.phase(4).integral;
guess.phase(6).integral = guess.phase(5).integral;
guess.phase(7).integral = guess.phase(6).integral;
guess.phase(8).integral = guess.phase(7).integral;
guess.phase(9).integral = guess.phase(7).integral;

%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 3;
mesh.colpointsmin = 2;
mesh.colpointsmax = 50;
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
setup.nlp.ipoptoptions.maxiterations = 1000;
% setup.derivatives.supplier           = 'sparseCD';
% setup.derivatives.derivativelevel    = 'second';
setup.derivatives.supplier           = 'sparseFD';
setup.derivatives.derivativelevel    = 'first';
setup.scales.method                  = 'automatic-bounds';
setup.method                         = 'RPM-Differentiation';
setup.scales.method                  = 'automatic-guessUpdate';

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

throttle = ones(length(alt),1);

throttle(time>output.result.solution.phase(1).time(end)) = 0.891;

throttle(time>output.result.solution.phase(2).time(end)) = 0.812;

throttle(time>output.result.solution.phase(3).time(end)) = .7333;

throttle(time>output.result.solution.phase(4).time(end)) = .6545;

throttle(time>output.result.solution.phase(5).time(end)) = .5757;

throttle(time>output.result.solution.phase(6).time(end)) = 0.496;

throttle(time>=output.result.solution.phase(7).time(end)) = 1; %after separation

throttle(time>=output.result.solution.phase(8).time(end)) = 0; %after separation


forward0 = [alt(1),gamma(1),v(1),zeta(1),lat(1),lon(1), mFuel(1)];


time_diff = [1 diff(time)];

time_forward = time;
Alpha_forward = Alpha;
throttle_forward = throttle;

time_forward(time_diff==0) = [];
Alpha_forward(time_diff==0) = [];
throttle_forward(time_diff==0) = [];

% stage = 1; %will be wrong for orbiter

stage_forward = ones(1,length(time_forward));
stage_forward(time_forward>=output.result.solution.phase(7).time(end)) = 2;


[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time_forward,Alpha_forward,f_t),ControlInterp(time_forward,throttle_forward,f_t),ControlInterp(time_forward,stage_forward,f_t)),0:time(end),forward0);
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
[altdot1,londot1,latdot1,gammadot1,vdot1(i),azidot1, q1(i), M(i), Fd(i), rho,L(i),Fueldt1,T(i),Isp1,Isp2,m(i),heating_rate(i),total_acceleration(i)] = SpaceLinerVehicleModel(time(i),phase_temp,throttle(i),auxdata,stage(i));

end

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
legend('Angle of Attack','Trajectory Angle')
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
plot(time,q1)
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


%%

    
    figure(231)
hold on

axesm('pcarree','Origin',[0 rad2deg(lon0) 0])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
plotm(rad2deg(lat),rad2deg(lon))
    
    
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
% 3: Check guess, is it unreasonable? Is it too close to the expected
% solution? Both can cause errors! Sometimes there is no real rhyme or
% reason to picking the correct guess, but a close bound to
% the expected solution has worked the most in my experience
% 4: Play with the no. of nodes, try both even and odd values
% 5: Play with scaling
% 6: Try all of the above in various combinations until it works!





