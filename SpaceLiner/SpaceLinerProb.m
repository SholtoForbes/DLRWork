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

Stage1.A = 93.25; %Reference Area in m� 
Stage1.mStruct = 215555.8;

Stage1.T_SL = 1961*9*1e3; %for mixture ratio 6 with all boosters active
Stage1.Isp_SL = 389;

Stage1.T_vac = 2206*9*1e3;
Stage1.Isp_vac = 437;

Stage2.A = 461; %Reference Area in m� 
Stage2.mStruct = 12800;

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
lat0 = -0.264;
lon0 = deg2rad(145);
aoaMin = 0;  aoaMax = 20*pi/180;
bankMin1 = -50*pi/180; bankMax1 =   50*pi/180;

% Primal Bounds
bounds.phase(1).state.lower = [0, lonMin, latMin, 0, -deg2rad(50), -2*pi, 0, aoaMin];
bounds.phase(1).state.upper = [100000, lonMax, latMax, 10000, deg2rad(89.9), 2*pi, 1.6038e+06, aoaMax];

% Initial States
bounds.phase(1).initialstate.lower = [1000,lon0, lat0, 100, deg2rad(89), deg2rad(60), 1.6038e+06, aoaMin] ;
bounds.phase(1).initialstate.upper = [1000,lon0, lat0, 100, deg2rad(89), deg2rad(60), 1.6038e+06, aoaMax];

% End States
bounds.phase(1).finalstate.lower = [80000, lonMin, latMin, 0, -deg2rad(70), -2*pi, 0, aoaMin];
bounds.phase(1).finalstate.upper = [80000, lonMax, latMax, 10000, deg2rad(70), 2*pi,1.6038e+06, aoaMax];

% Control Bounds
bounds.phase(1).control.lower = [deg2rad(-.1)];
bounds.phase(1).control.upper = [deg2rad(.1)];

% Time Bounds
bounds.phase(1).initialtime.lower = 0;
bounds.phase(1).initialtime.upper = 0;
bounds.phase(1).finaltime.lower = 100;
bounds.phase(1).finaltime.upper = 5000;

%% Define Path Constraints
% Path bounds, defined in Continuous function.
% These limit the dynamic pressure.

bounds.phase(1).path.lower = [0];
bounds.phase(1).path.upper = [60000];


%%  Guess =================================================================
% Set the initial guess. This can have a significant effect on the final
% solution, even for a well defined problem. 
guess.phase(1).state(:,1)   = [2000;70000];
guess.phase(1).state(:,2)   = [0;0];
guess.phase(1).state(:,3)   = [-0.269;-0.13];
guess.phase(1).state(:,4)   = [100,10000];
guess.phase(1).state(:,5)   = [deg2rad(80),0];
guess.phase(1).state(:,6)   = [deg2rad(60),deg2rad(60)];
guess.phase(1).state(:,7) 	= [1.6038e+06, 100];
guess.phase(1).state(:,8)   = [15*pi/180; 15*pi/180];
% guess.phase(1).state(:,9)   = [deg2rad(0);deg2rad(0)];


guess.phase(1).control      = [[0;0]];
guess.phase(1).time          = [0;650];

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

eta = output.result.solution.phase(1).state(:,9).';





time = output.result.solution.phase(1).time.';


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


% =========================================================================

%% Third Stage
% Optimise third stage trajectory from end point

global phi



if auxdata.PayloadGrid(alt(end)+10,gamma(end),v(end)) - auxdata.PayloadGrid(alt(end),gamma(end),v(end)) < 0
    disp('Check Third Stage Payload Matrix, Found Maxima')
end

if gamma(end) > max(ThirdStageData(:,4)) || gamma(end) < min(ThirdStageData(:,4))
    disp('Third Stage Matrix Extrapolating for Trajectory Angle')
end

if alt(end) > max(ThirdStageData(:,3)) || alt(end) < min(ThirdStageData(:,3))
    disp('Third Stage Matrix Extrapolating for Altitude')
end

if v(end) > max(ThirdStageData(:,5)) || v(end) < min(ThirdStageData(:,5))
    disp('Third Stage Matrix Extrapolating for Velocity')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          OUTPUT             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes = length(alt)

% eq = Engine.eq;
% Thrust = Engine.Thrust;
% Fueldt = Engine.Fueldt;

% Fd = Vehicle.Fd;
% Alpha = Vehicle.Alpha;
% lift = Vehicle.lift;
% flapdeflection = Vehicle.flapdeflection;
% 
% Thrust = Thrust./cos(deg2rad(Alpha)); % change thrust to account for total thrust, including portion that contributes to lift
% 
% dt = time(2:end)-time(1:end-1); % Time change between each node pt
% FuelUsed = zeros(1,nodes-1);
% FuelUsed(1) = dt(1)*Fueldt(1);
% for i = 2:nodes-1
%     FuelUsed(i) = dt(i).*Fueldt(i) + FuelUsed(i-1);
% end



[~,~,~,~,~,~, q1, M1, Fd1, rho1,L1,Fueldt1,T1,Isp1,q11,flapdeflection] = VehicleModelCombined(gamma, alt, v,auxdata,zeta,lat,lon,Alpha,eta,1, mFuel,mFuel(1),mFuel(end), 1);
[~,~,~,~,~,~, q2, M2, Fd2, rho2,L2,Fueldt2,T2,Isp2,q12,flapdeflection2] = VehicleModelCombined(gamma2, alt2, v2,auxdata,zeta2,lat2,lon2,Alpha2,eta2,throttle2, mFuel2,0,0, 0);

throttle2(M2<5.0) = 0; % remove nonsense throttle points

% figure out horizontal motion
H(1) = 0;
for i = 1:nodes-1
H(i+1) = v(i)*(time(i+1) - time(i))*cos(gamma(i)) + H(i);
end

% Separation_LD = lift(end)/Fd(end)

figure(2010)

subplot(5,5,[1,10])
hold on
plot(H, alt)
% plot(H(algorithm.nodes(1)), V(algorithm.nodes(1)), '+', 'MarkerSize', 10, 'MarkerEdgeColor','r')
title('Trajectory (m)')

dim = [.7 .52 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass), ' kg'],['Second Stage Fuel Used: ' num2str(1562 - mFuel(end)) ' kg']},'FitBoxToText','on');  


subplot(5,5,11)
hold on
plot(time, v)

title('Velocity (m/s)')


subplot(5,5,12)
plot(time, M1)
title('Mach no')

subplot(5,5,13)
plot(time, q1)
title('Dynamic Pressure (pa)')

subplot(5,5,14)
hold on
plot(time, rad2deg(gamma))

title('Trajectory Angle (Deg)')



subplot(5,5,15)
plot(time, Fd1)
title('Drag Force')

subplot(5,5,16)
hold on
plot(time, mFuel + 8755.1 - 994)
title('Vehicle Mass (kg)')



subplot(5,5,17)
plot(time, T1)
title('Thrust (N)')

% Isp1 = T1./Fueldt1./9.81;
IspNet1 = (T1-Fd1)./Fueldt1./9.81;

subplot(5,5,18)
plot(time, Isp1)
title('Isp')

subplot(5,5,19)
plot(time, IspNet1)
title('Net Isp')

subplot(5,5,20)
plot(time, flapdeflection)
title('Flap Deflection (deg)')

subplot(5,5,21)
plot(time, rad2deg(Alpha))
title('Angle of Attack (deg)')

subplot(5,5,22)
plot(time, rad2deg(eta))
title('Bank Angle (deg)')

subplot(5,5,23)
plot(time, q11)
title('Dynamic pressure after shock')

% subplot(5,5,22);
% plot(time, dual.dynamics);
% title('costates')
% xlabel('time');
% ylabel('costates');
% legend('\lambda_1', '\lambda_2', '\lambda_3');

% subplot(5,5,23)
% Hamiltonian = dual.Hamiltonian(1,:);
% plot(time,Hamiltonian);
% title('Hamiltonian')

% subplot(5,5,24)
% hold on
% plot(time, rad2deg(gammadot))
% title('Trajectory Angle Change Rate (Deg/s)')
% 
% subplot(5,5,25)
% hold on
% plot(time, rad2deg(omegadot))
% title('Omegadot Control (Deg/s2)')


dim = [.8 .0 .2 .2];
annotation('textbox',dim,'string',{['Third Stage Thrust: ', num2str(50), ' kN'],['Third Stage Starting Mass: ' num2str(2850) ' kg'],['Third Stage Isp: ' num2str(350) ' s']},'FitBoxToText','on');  

figure(202)
sp1 = subplot(2,6,[1,6]);
ax1 = gca; % current axes
hold on
plot(H/1000, alt/1000,'Color','k')

title('Trajectory')
xlabel('Earth Normal Distance Flown (km)')
ylabel('Vertical Position (km)')

for i = 1:floor(time(end)/30)
    [j,k] = min(abs(time-30*i));
    str = strcat(num2str(round(time(k))), 's');
    text(H(k)/1000,alt(k)/1000,str,'VerticalAlignment','top', 'FontSize', 10);
    
    plot(H(k)/1000, alt(k)/1000, '+', 'MarkerSize', 10, 'MarkerEdgeColor','k')
end

plot(H(end)/1000, alt(end)/1000, 'o', 'MarkerSize', 10, 'MarkerEdgeColor','k')

text(H(end)/1000,alt(end)/1000,'Third Stage Transition Point','VerticalAlignment','top', 'FontSize', 10);

dim = [.65 .45 .2 .2];
annotation('textbox',dim,'string',{['Payload Mass: ', num2str(ThirdStagePayloadMass,4), ' kg'],['Second Stage Fuel Used: ' num2str(mFuel(1) - mFuel(end)) ' kg']},'FitBoxToText','on');  

thirdstageexample_H = [0+H(end) (H(end)-H(end - 1))+H(end) 20*(H(end)-H(end - 1))+H(end) 40*(H(end)-H(end - 1))+H(end) 60*(H(end)-H(end - 1))+H(end) 80*(H(end)-H(end - 1))+H(end)]/1000; %makes a small sample portion of an arbitrary third stage trajectory for example
thirdstageexample_V = [0+alt(end) (alt(end)-alt(end - 1))+alt(end) 20*((alt(end)-alt(end -1)))+alt(end) 40*((alt(end)-alt(end -1)))+alt(end) 60*((alt(end)-alt(end -1)))+alt(end) 80*((alt(end)-alt(end -1)))+alt(end)]/1000;
plot(thirdstageexample_H, thirdstageexample_V, 'LineStyle', '--','Color','k');

hold on
sp2 = subplot(2,6,[7,9]);
xlabel('time (s)')

hold on
ax2 = gca; % current axes
xlim([min(time) max(time)]);

line(time, rad2deg(gamma),'Parent',ax2,'Color','k', 'LineStyle','-')

line(time, M1,'Parent',ax2,'Color','k', 'LineStyle','--')

line(time, v./(10^3),'Parent',ax2,'Color','k', 'LineStyle','-.')

line(time, q1./(10^4),'Parent',ax2,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)

% line(time, heating_rate./(10^5),'Parent',ax1,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% 
% line(time, Q./(10^7),'Parent',ax1,'Color','k', 'LineStyle','-', 'lineWidth', 2.0)

% legend(ax1,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)',  'Q (Mj x 10)')
h = legend(ax2,  'Trajectory Angle (degrees)', 'Mach no', 'Velocity (m/s x 10^3)', 'Dynamic Pressure (Pa x 10^4)');
rect1 = [0.12, 0.35, .25, .25];
set(h, 'Position', rect1)


sp3 = subplot(2,6,[10,12]);
xlabel('time (s)')
ax3 = gca;
xlim([min(time) max(time)]);
line(time, [rad2deg(Alpha(1:end-1)) rad2deg(Alpha(end-1))],'Parent',ax3,'Color','k', 'LineStyle','-')

line(time, flapdeflection,'Parent',ax3,'Color','k', 'LineStyle','--')


% line(time, mfuel./(10^2),'Parent',ax2,'Color','k', 'LineStyle','-.')
% line(time, eq.*10,'Parent',ax3,'Color','k', 'LineStyle','-.')

line(time, IspNet1./(10^2),'Parent',ax3,'Color','k', 'LineStyle',':', 'lineWidth', 2.0)
% 
% g = legend(ax2, 'AoA (degrees)','Flap Deflection (degrees)', 'Fuel Mass (kg x 10^2)', 'Net Isp (s x 10^2)');
g = legend(ax3, 'AoA (degrees)','Flap Deflection (degrees)', 'Equivalence Ratio x 10', 'Net Isp (s x 10^2)');

rect2 = [0.52, 0.35, .25, .25];
set(g, 'Position', rect2)

saveas(figure(202),[sprintf('../ArchivedResults/%s',Timestamp),filesep,'SecondStage.fig']);




%% PLOT RETURN
addpath('..\SecondStageReturn\addaxis')

figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% subplot(3,1,1)
hold on


 plot(time2,alt2/1000,'-','color','k', 'linewidth', 1.);
% ylim([-30,40]);
ylabel('altitude(km)');
xlabel('time (s)');
addaxis(time2,v2/1000,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Velocity (km/s)');

% addaxis(time,fpa,':','color','k', 'linewidth', 1.);
% addaxislabel(3,'Trajectory Angle (deg)');

addaxis(time2,zeta2,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Heading Angle (Deg)');


legend(  'Altitude', 'Velocity', 'Heading Angle', 'location', 'best');

figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
% subplot(3,1,2)
hold on
plot(time2,rad2deg(Alpha2),'-','color','k', 'linewidth', 1.);
ylabel('Angle of Attack (deg)');
xlabel('time (s)')
throttle(M2<5.0)=0;
addaxis(time2,throttle2*100,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Throttle (%)');

% addaxis(time,mfuel,'-.','color','k', 'linewidth', 1.);
% addaxislabel(3,'Fuel Mass (kg)');



addaxis(time2,rad2deg(eta2),':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Bank Angle (Deg)');

addaxis(time2,flapdeflection2,'-.','color','k', 'linewidth', 1.2);
addaxislabel(4,'Flap Deflection (Deg)');
legend(  'Angle of Attack', 'Throttle' , 'Bank Angle','FlapDeflection');


figure('units','normalized','outerposition',[0.1 0.1 .7 .5])
hold on
% subplot(3,1,3)
plot(time2,M2,'-','color','k', 'linewidth', 1.);
ylabel('Mach no.')
xlabel('time (s)')

addaxis(time2,Isp2,'--','color','k', 'linewidth', 1.);
addaxislabel(2,'Specific Impulse (s)');

addaxis(time2,q2,':','color','k', 'linewidth', 1.2);
addaxislabel(3,'Dynamic Pressure (kPa)');

% addaxis(time,L./D,':','color','k', 'linewidth', 1.);
% addaxislabel(4,'L/D');

legend(  'Mach no.', 'Isp (potential)', 'q' );


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FORWARD SIMULATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% This is a full forward simulation, using the angle of attack and flap
% deflection at each node.

% Note, because the nodes are spaced widely, small interpolation
% differences result in the forward simulation being slightly different
% than the actual. This is mostly a check to see if they are close. 


forward0 = [alt(1),gamma(1),v(1),zeta(1),lat(1),lon(1), mFuel(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelAscent_forward(f_t, f_y,auxdata,ControlInterp(time,Alpha,f_t),ControlInterp(time,eta,f_t),1,mFuel(1),mFuel(end)),time(1:end),forward0);

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



% Return Forward
forward0 = [alt2(1),gamma2(1),v2(1),zeta2(1),lat2(1),lon2(1), mFuel2(1)];

% [f_t, f_y] = ode45(@(f_t,f_y) ForwardSim(f_y,AlphaInterp(t,Alpha,f_t),communicator,communicator_trim,SPARTAN_SCALE,Atmosphere,const,scattered),t,forward0);
[f_t, f_y] = ode45(@(f_t,f_y) VehicleModelReturn_forward(f_t, f_y,auxdata,ControlInterp(time2,Alpha2,f_t),ControlInterp(time2,eta2,f_t),ControlInterp(time2,throttle2,f_t)),time2(1:end),forward0);

% altitude  = (output.result.solution.phase(1).state(:,1)-auxdata.Re);
figure(213)
subplot(7,1,1)
hold on
plot(f_t(1:end),f_y(:,1));
plot(time2,alt2);

% gamma  = output.result.solution.phase.state(:,5);

subplot(7,1,2)
hold on
plot(f_t(1:end),f_y(:,2));
plot(time2,gamma2);

% latitude  = output.result.solution.phase.state(:,3);
subplot(7,1,3:5)
hold on
plot(f_y(:,6),f_y(:,5));
plot(lon2,lat2);

subplot(7,1,6)
hold on
plot(f_t(1:end),f_y(:,7));
plot(time2,mFuel2);

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





