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

auxdata.interp.Stage1.Cl_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cl_Grid','spline');
auxdata.interp.Stage1.Cd_spline = griddedInterpolant(Stage1_MGrid',Stage1_aoaGrid',Stage1_Cd_Grid','spline');

auxdata.interp.Stage2.Cl_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cl_Grid','spline');
auxdata.interp.Stage2.Cd_spline = griddedInterpolant(Stage2_MGrid',Stage2_aoaGrid',Stage2_Cd_Grid','spline');


%% Import Bounds %%========================================================
lonMin = -pi;         lonMax = -lonMin;
latMin = -70*pi/180;  latMax = -latMin;
lat0 = deg2rad(45);
lon0 = deg2rad(145);
aoaMin = 0;  aoaMax = 20*pi/180;
% bankMin1 = -50*pi/180; bankMax1 =   50*pi/180;

% Primal Bounds
bounds.phase(1).state.lower = [0, lonMin, latMin, 0, -deg2rad(50), -pi, 0, aoaMin];
bounds.phase(1).state.upper = [200000, lonMax, latMax, 15000, deg2rad(87), pi, 1.6038e+06, aoaMax];

% Initial States
% bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 95, deg2rad(80), deg2rad(70), 1.5038e+06, aoaMin] ;
% bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 105, deg2rad(89), deg2rad(80), 1.6038e+06, aoaMax];

bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 95, deg2rad(80), deg2rad(70), 1.28e+06, aoaMin] ;
bounds.phase(1).initialstate.upper = [1200,lon0, lat0, 105, deg2rad(87), deg2rad(80), 1.32e+06, aoaMax];

% End States

% bounds.phase(1).finalstate.lower = [70000, lonMin, latMin, 0, -deg2rad(60), -2*pi, 0, aoaMin];
% bounds.phase(1).finalstate.upper = [70000, lonMax, latMax, 10000, deg2rad(60), 2*pi,1.6038e+06, aoaMax];
% 
bounds.phase(1).finalstate.lower = [70000, lonMin, latMin, 0, 0, -2*pi, 0, aoaMin];
bounds.phase(1).finalstate.upper = [70000, lonMax, latMax, 10000, 0, 2*pi,1.6038e+06, aoaMax];


% Control Bounds
bounds.phase(1).control.lower = [deg2rad(-.5)];
bounds.phase(1).control.upper = [deg2rad(.5)];

% Time Bounds
bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower = 100;
bounds.phase(1).finaltime.upper = 5000;

%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.

bounds.phase(1).path.lower = [0];
bounds.phase(1).path.upper = [40000];

%% Bound integral if necessary
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 1e9;

%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(1).state(:,1)   = [2000;70000];
guess.phase(1).state(:,2)   = [2.5;2.7];
guess.phase(1).state(:,3)   = [0.7;1];
guess.phase(1).state(:,4)   = [100,10000];
guess.phase(1).state(:,5)   = [deg2rad(80),0];
guess.phase(1).state(:,6)   = [deg2rad(70),deg2rad(70)];
guess.phase(1).state(:,7) 	= [1.6038e+06, 1e5];
guess.phase(1).state(:,8)   = [1*pi/180; 1*pi/180];
% guess.phase(1).state(:,9)   = [deg2rad(0);deg2rad(0)];


guess.phase(1).control      = [[0;0]];
guess.phase(1).time          = [0;450];


guess.phase.integral = 0;
% Tire stages together
% bounds.eventgroup(1).lower = [zeros(1,10)];
% bounds.eventgroup(1).upper = [zeros(1,10)]; 

%%
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
% mesh.method       = 'hp-LiuRao-Legendre';
mesh.maxiterations = 10;
mesh.colpointsmin = 3;
mesh.colpointsmax = 100;
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

mFuel = output.result.solution.phase(1).state(:,7).'; 

Alpha = output.result.solution.phase(1).state(:,8).';

% eta = output.result.solution.phase(1).state(:,9).';

states.alt     = output.result.solution.phase(1).state(:,1);
states.lon     = output.result.solution.phase(1).state(:,2);
states.lat     = output.result.solution.phase(1).state(:,3);
states.v       = output.result.solution.phase(1).state(:,4);
states.gamma   = output.result.solution.phase(1).state(:,5);
states.zeta    = output.result.solution.phase(1).state(:,6);
states.mFuel   = output.result.solution.phase(1).state(:,7);



time = output.result.solution.phase(1).time.';

controls.Alpha = Alpha';

figure(201)
subplot(9,1,1)
hold on
plot(time,alt)
subplot(9,1,2)
hold on
plot(time,v)
subplot(9,1,3)
hold on
plot(time,lon)
subplot(9,1,4)
hold on
plot(time,lat)
subplot(9,1,5)
hold on
plot(time,v)
subplot(9,1,6)
hold on
plot(time,gamma)
subplot(9,1,7)
hold on
plot(time,ones(1,length(time)))



figure(230)
hold on
plot3(lon,lat,alt)


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 

throttle = ones(length(alt),1);

throttle(time>time(end)/20*5) = 0.891;

throttle(time>time(end)/20*6) = 0.812;

throttle(time>time(end)/20*7) = .7333;

throttle(time>time(end)/20*8) = .6545;

throttle(time>time(end)/20*9) = .5757;

throttle(time>time(end)/20*10) = 0.496;

throttle(time>=time(end)*0.6) = 1; %after separation

forward0 = [alt(1),gamma(1),v(1),zeta(1),lat(1),lon(1), mFuel(1)];

[f_t, f_y] = ode45(@(f_t,f_y) VehicleModel_forward(f_t, f_y,auxdata,ControlInterp(time,Alpha,f_t),ControlInterp(time,throttle,f_t),time(end)),time(1:end),forward0);

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



[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1, M, Fd, rho,L,Fueldt1,T,Isp1,Isp2,m] = SpaceLinerVehicleModel(time,states,controls,throttle,auxdata,time(end));

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





