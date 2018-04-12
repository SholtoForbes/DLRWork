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

Stage1.A = 461;% %Reference Area in m² 


Stage1.mStruct = 205245.8;


auxdata.Stage1 = Stage1;


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

Stage1_MList = [ 0.10  0.35  0.60  0.85  1.10  1.35  1.60  1.85  2.10  2.35  2.60  2.85  3.10 4.000  6.000  8.000 10.000 12.000 14.000];
Stage1_aoaList = [-10.00  -8.00  -6.00  -4.00  -2.00   0.00   2.00   4.00   6.00   8.00  10.00  12.00  14.00  16.00  18.00  20.00  22.00  24.00  26.00  28.00  30.00  32.00  34.00  36.00  38.00  40.00  42.00];
Stage1_Cl = importdata('aero_Booster_Lift');
Stage1_Cd = importdata('aero_Booster_Drag');

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
        Stage1_Cl_Grid(j,i) = Stage1_Cl_scattered(Stage1_MGrid(j,i),Stage1_aoaGrid(j,i));
        Stage1_Cd_Grid(j,i) = Stage1_Cd_scattered(Stage1_MGrid(j,i),Stage1_aoaGrid(j,i));
    end
end

%% Create Gridded Interpolants

auxdata.interp.Stage1.Cl_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cl_Grid','linear');
auxdata.interp.Stage1.Cd_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cd_Grid','linear');


%% Get land interpolants to determine if over land
load LandSpline

auxdata.LandSpline = LandSpline;


%% Get population interpolation
load PopInterp

auxdata.PopInterp = PopInterp;

%% Set Bounds %%========================================================



%%

altMin = 1;
altMax = 100000;


vMin = 1;
vMax = 10000;

gammaMin = -deg2rad(80);
% gammaMax = deg2rad(80);
gammaMax = deg2rad(80);

zetaMin = -3*pi;
% zetaMin = 0;
zetaMax = 3*pi;

lonMin = -2*pi;         
lonMax = 2*pi;

latMin = -pi/2+0.0000001;  
% latMin = 0;  %for japan-germany actual

latMax = pi/2-0.0000001;
% latMax = pi;

aoaMin = deg2rad(-10);  
aoaMax = deg2rad(40);

bankMin = deg2rad(-80); 
bankMax =   deg2rad(80);


%% Initial Conditions
mission = 1;

if mission == 4
%aus-germany
alt0 = 41002;
v0 = 3726.4;
gamma0 = -0.0103;
lat0 = -0.3431; % 
lon0 = 0.0082+deg2rad(150.809994); % 
zeta0 = 1.4247; %
end

if mission == 1
% Japan-Germany
alt0 = 48484;
v0 = 3955.8;
gamma0 = 0.0884;
lat0 = 0.7186; % 
lon0 = 0.0087+deg2rad(137.276039); % 
zeta0 = 1.4895; %
end

if mission == 2
%Germany-Japan
alt0 = 42720;
v0 = 3648.3;
gamma0 = 0.0166;
lat0 = 1.0038; % 
lon0 = -0.0419+deg2rad(7.763275); % 
zeta0 = 1.9264; %
end

if mission ==3
%Florida-Aus
alt0 = 43196;
v0 = 3862.4;
gamma0 = 0.0410;
lat0 = 0.5071; % 
lon0 = -0.0739+deg2rad(-83.159933); % 
zeta0 = 3.2529; %
end

if mission == 5
%Aus-Florida
alt0 = 44199;
v0 = 3720.3;
gamma0 = 0.0513;
lat0 = -0.3875; % 
lon0 =0.07+deg2rad(150.809994); % 
zeta0 = 0.3070; %
end



auxdata.lon0 = lon0;

altFMin = 5000;
altFMax = 10000;

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

% latF = deg2rad(37.235705); %south korea
% lonF = deg2rad(129.356220);

% latF = deg2rad(38.745095); %north korea
% lonF = deg2rad(128.272375);


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

bounds.phase(1).initialstate.lower = [alt0,0, lat0, v0, gamma0, zeta0, aoaMin, 0] ;
bounds.phase(1).initialstate.upper = [alt0,0, lat0, v0, gamma0, zeta0, aoaMax, 0];


% End States
%Starting East heading east
% For aus-germany, japan-germany
% if mission == 3 || mission == 2
% zetaF = zeta0+pi;
% elseif mission == 1 
%    zetaF = zeta0-pi; 
% end

bounds.phase(1).finalstate.lower = [altFMin, lonMin, latMin, 230, deg2rad(-20), zetaMin, aoaMin, 0];
bounds.phase(1).finalstate.upper = [altFMax, lonMax, latMax, 235, deg2rad(20), zetaMax, aoaMax, 0];


% Control Bounds
bounds.phase(1).control.lower = [deg2rad(-1), deg2rad(-3)];
bounds.phase(1).control.upper = [deg2rad(1), deg2rad(3)];

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

bounds.phase(1).path.lower = [0, 0, -10];
bounds.phase(1).path.upper = [40000, 1.3e6, 10];

% bounds.phase(1).path.lower = [0, 0, -1.5]; 
% bounds.phase(1).path.upper = [30000, 1.3e6, 1.5];

% bounds.phase(1).path.lower = [0, 0, -1]; 
% bounds.phase(1).path.upper = [30000, 1.5e6, 1];

% bounds.phase(1).path.lower = [0, 0, -1]; 
% bounds.phase(1).path.upper = [7000, 1.3e6, 1];

%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(1).state(:,1)   = [alt0;alt0];
guess.phase(1).state(:,2)   = [0;0]; %East
% guess.phase(1).state(:,2)   = [0;lonF-lon0-2*pi]; %West

% guess.phase(1).state(:,2)   = [0;lonF-lon0];
guess.phase(1).state(:,3)   = [lat0;lat0];
% guess.phase(1).state(:,3)   = [lat0;latF-0.5];
guess.phase(1).state(:,4)   = [v0;vMin];

% This is gamma
guess.phase(1).state(:,5)   = [0;0]; %


guess.phase(1).state(:,6)   = [zeta0;zeta0]; % Aus-Germany

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
mesh.maxiterations = 4;
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
if mission == 1
setup.nlp.ipoptoptions.maxiterations = 900;

elseif mission ==2
  setup.nlp.ipoptoptions.maxiterations = 700;  

elseif mission ==3
  setup.nlp.ipoptoptions.maxiterations = 700; 
  
  elseif mission ==4
  setup.nlp.ipoptoptions.maxiterations = 800;  
    
  elseif mission ==5
  setup.nlp.ipoptoptions.maxiterations = 600;  
    
  
end

setup.derivatives.supplier           = 'sparseCD';
% setup.derivatives.derivativelevel    = 'second';
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

[altdot,londot,latdot,gammadot,a,zetadot, q, M, D, rho,L,Fueldt,heating_rate,az] = SpaceLinerVehicleModel(states,controls,auxdata);

total_acceleration = sqrt(a.^2 + (states.v.*gammadot).^2 + (states.v.*zetadot).^2)/9.81;
% total_acceleration = a/9.81;

pop = auxdata.PopInterp(rad2deg(states.lon+auxdata.lon0),rad2deg(states.lat));


figure(201)
subplot(5,2,1)
hold on
plot(time,alt)
xlabel('time (s)')
ylabel('Altitude (km)')

subplot(5,2,2)
hold on
plot(v,alt/1000)
xlabel('Velocity (m/s)')
ylabel('Altitude (km)')


subplot(5,2,3)
hold on
plot(time,rad2deg(Alpha))
plot(time,rad2deg(gamma))
plot(time,rad2deg(eta))
xlabel('time (s)')
ylabel('(deg)');
legend('Angle of Attack' ,'Trajectory Angle', 'Bank Angle')

subplot(5,2,4)
hold on
plot(time,D/1000)
plot(time,L/1000)
legend('Drag', 'Lift')
xlabel('time (s)')
ylabel('(kN)');

subplot(5,2,5)
hold on
plot(time,abs(total_acceleration))
xlabel('time (s)')
ylabel('Acceleration (g)')

subplot(5,2,6)
hold on
plot(time,heating_rate/1e6)
xlabel('time (s)')
ylabel('Heating Rate (MW/m^2')


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





figure(204)
subplot(2,2,1)
hold on
plot(time,alt/1000,'k')
xlabel('time (s)')
ylabel('altitude (km)')

subplot(2,2,2)
hold on
plot(v,alt/1000,'k')
xlabel('velocity (m/s)')
ylabel('altitude (km)')


subplot(2,2,3)
hold on
plot(time,rad2deg(Alpha),'k')
plot(time,rad2deg(gamma))
plot(time,rad2deg(eta))
legend('Angle of Attack','Trajectory Angle', 'Bank Angle')
xlabel('Time (s)')
ylabel('(deg)')

subplot(2,2,4)
hold on

if mission == 2
xlim([rad2deg(lon0)-15,rad2deg(lon0)+3])
ylim([rad2deg(lat0)-3,rad2deg(lat0)+7])

elseif mission == 4 
xlim([rad2deg(lon0)-3,rad2deg(lon0)+7])
ylim([rad2deg(lat0)-3,rad2deg(lat0)+7])

elseif   mission == 5 
xlim([rad2deg(lon0)-3,rad2deg(lon0)+10])
ylim([rad2deg(lat0)-3,rad2deg(lat0)+7])

elseif mission ==1 
xlim([rad2deg(lon0)-3,rad2deg(lon0)+7])
ylim([rad2deg(lat0)-1,rad2deg(lat0)+9])

elseif mission == 3
xlim([rad2deg(lon0)-10,rad2deg(lon0)+3])
ylim([rad2deg(lat0)-4,rad2deg(lat0)+5])
end

% axesm('pcarree','Origin',[0 rad2deg(lon0)+60 0])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
plot(rad2deg(lon)+rad2deg(lon0),rad2deg(lat),'k','LineWidth',2)


    xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
    

    
figure(205)

subplot(2,2,1)
hold on
plot(time,q,'k')
xlabel('time (s)')
ylabel('Dynamic Pressure (kPa)')


subplot(2,2,4)
hold on

plot(time,az/9.81,'k')

xlabel('time (s)')
ylabel('acceleration (g)')


subplot(2,2,3)
hold on
plot(M,alt/1000,'k')
xlabel('Mach no.')
ylabel('Altitude (km)')

subplot(2,2,2)
hold on
plot(time,heating_rate/1e6,'k')

xlabel('time (s)')
ylabel('Heating Rate (MW/m2)')



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
    
    
%% save file
contents = {'time', 'v', 'gamma', 'zeta', 'alt', 'lon', 'lat', 'q', 'heating_rate', 'L', 'D', 'M', 'Alpha' ,'eta'};
units = {'s', 'm/s', 'rad', 'rad', 'm', 'rad', 'rad', 'pa', 'W/m^2', 'N', 'N', 'Mach', 'rad', 'rad'};

result = [time' v' gamma' zeta' alt' lon' lat' total_acceleration q heating_rate L D M Alpha' eta' ];
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
