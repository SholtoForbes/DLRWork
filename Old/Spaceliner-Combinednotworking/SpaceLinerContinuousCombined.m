function phaseout = SpaceLinerContinuousCombined(input)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Dynamics


auxdata = input.auxdata;


% throttle = ones(length(states.alt),1);

% throttle(time>time(end)/20*5) = 0.891; %using throttle to turn engines off
% 
% throttle(time>time(end)/20*6) = 0.812;
% 
% throttle(time>time(end)/20*7) = .7333;
% 
% throttle(time>time(end)/20*8) = .6545;
% 
% throttle(time>time(end)/20*9) = .5757;
% 
% throttle(time>time(end)/20*10) = 0.496;
% 
% throttle(time>=time(end)*0.6) = 1; %after separation

%% 1
Alphadot1  = input.phase(1).control(:,1)';
etadot1  = input.phase(1).control(:,2)';
time1 = input.phase(1).time;
throttle1 = 1;
[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1, M1, Fd1, rho1,L1,Fueldt1,T1,Isp1,Isp2,m1,heating_rate1,total_acceleration1] = SpaceLinerVehicleModel(time1,input.phase(1),throttle1,auxdata,1);
phaseout(1).dynamics  = [altdot1', londot1', latdot1', vdot1', gammadot1', azidot1', -Fueldt1, Alphadot1', etadot1'];
phaseout(1).path = [q1,total_acceleration1'];
phaseout(1).integrand = heating_rate1;

%% 2
Alphadot2  = input.phase(2).control(:,1)';
etadot2  = input.phase(2).control(:,2)';
time2 = input.phase(2).time;
throttle2 = 0.891;
[altdot2,londot2,latdot2,gammadot2,vdot2,azidot2, q2, M2, Fd2, rho2,L2,Fueldt2,T2,Isp2,Isp2,m2,heating_rate2,total_acceleration2] = SpaceLinerVehicleModel(time2,input.phase(2),throttle2,auxdata,1);
phaseout(2).dynamics  = [altdot2', londot2', latdot2', vdot2', gammadot2', azidot2', -Fueldt2, Alphadot2', etadot2'];
phaseout(2).path = [q2,total_acceleration2'];
phaseout(2).integrand = heating_rate2;

%% 3
Alphadot3  = input.phase(3).control(:,1)';
etadot3  = input.phase(3).control(:,2)';
time3 = input.phase(3).time;
throttle3 = 0.812;
[altdot3,londot3,latdot3,gammadot3,vdot3,azidot3, q3, M3, Fd3, rho3,L3,Fueldt3,T3,Isp3,Isp3,m3,heating_rate3,total_acceleration3] = SpaceLinerVehicleModel(time3,input.phase(3),throttle3,auxdata,1);
phaseout(3).dynamics  = [altdot3', londot3', latdot3', vdot3', gammadot3', azidot3', -Fueldt3, Alphadot3', etadot3'];
phaseout(3).path = [q3,total_acceleration3'];
phaseout(3).integrand = heating_rate3;

%% 4
Alphadot4  = input.phase(4).control(:,1)';
etadot4  = input.phase(4).control(:,2)';
time4 = input.phase(4).time;
throttle4 = .7333;
[altdot4,londot4,latdot4,gammadot4,vdot4,azidot4, q4, M4, Fd4, rho4,L4,Fueldt4,T4,Isp4,Isp4,m4,heating_rate4,total_acceleration4] = SpaceLinerVehicleModel(time4,input.phase(4),throttle4,auxdata,1);
phaseout(4).dynamics  = [altdot4', londot4', latdot4', vdot4', gammadot4', azidot4', -Fueldt4, Alphadot4', etadot4'];
phaseout(4).path = [q4,total_acceleration4'];
phaseout(4).integrand = heating_rate4;

%% 5
Alphadot5  = input.phase(5).control(:,1)';
etadot5  = input.phase(5).control(:,2)';
time5 = input.phase(5).time;
throttle5 = .6545;
[altdot5,londot5,latdot5,gammadot5,vdot5,azidot5, q5, M5, Fd5, rho5,L5,Fueldt5,T5,Isp5,Isp5,m5,heating_rate5,total_acceleration5] = SpaceLinerVehicleModel(time5,input.phase(5),throttle5,auxdata,1);
phaseout(5).dynamics  = [altdot5', londot5', latdot5', vdot5', gammadot5', azidot5', -Fueldt5, Alphadot5', etadot5'];
phaseout(5).path = [q5,total_acceleration5'];
phaseout(5).integrand = heating_rate5;

%% 6
Alphadot6  = input.phase(6).control(:,1)';
etadot6  = input.phase(6).control(:,2)';
time6 = input.phase(6).time;
throttle6 = .5757;
[altdot6,londot6,latdot6,gammadot6,vdot6,azidot6, q6, M6, Fd6, rho6,L6,Fueldt6,T6,Isp6,Isp6,m6,heating_rate6,total_acceleration6] = SpaceLinerVehicleModel(time6,input.phase(6),throttle6,auxdata,1);
phaseout(6).dynamics  = [altdot6', londot6', latdot6', vdot6', gammadot6', azidot6', -Fueldt6, Alphadot6', etadot6'];
phaseout(6).path = [q6,total_acceleration6'];
phaseout(6).integrand = heating_rate6;

%% 7
Alphadot7  = input.phase(7).control(:,1)';
etadot7  = input.phase(7).control(:,2)';
time7 = input.phase(7).time;
throttle7 = 0.496;
[altdot7,londot7,latdot7,gammadot7,vdot7,azidot7, q7, M7, Fd7, rho7,L7,Fueldt7,T7,Isp7,Isp7,m7,heating_rate7,total_acceleration7] = SpaceLinerVehicleModel(time7,input.phase(7),throttle7,auxdata,1);
phaseout(7).dynamics  = [altdot7', londot7', latdot7', vdot7', gammadot7', azidot7', -Fueldt7, Alphadot7', etadot7'];
phaseout(7).path = [q7,total_acceleration7'];
phaseout(7).integrand = heating_rate7;

%% 8
Alphadot8  = input.phase(8).control(:,1)';
etadot8  = input.phase(8).control(:,2)';
time8 = input.phase(8).time;
throttle8 = 1;
[altdot8,londot8,latdot8,gammadot8,vdot8,azidot8, q8, M8, Fd8, rho8,L8,Fueldt8,T8,Isp8,Isp8,m8,heating_rate8,total_acceleration8] = SpaceLinerVehicleModel(time8,input.phase(8),throttle8,auxdata,2);
% length(altdot8), length(londot8), length(latdot8), length(vdot8), length(gammadot8), length(azidot8), length(Fueldt8),length( Alphadot8), length(etadot8)
phaseout(8).dynamics  = [altdot8', londot8', latdot8', vdot8', gammadot8', azidot8', -Fueldt8, Alphadot8', etadot8'];
phaseout(8).path = [q8,total_acceleration8'];
phaseout(8).integrand = heating_rate8;


%% Descent
% mFuel8 = input.phase(8).state(:,7);
% auxdata.mFuel_descent = mFuel8(end);


% 2D Dynamics
states.alt9     = input.phase(9).state(:,1);
states.lon9     = input.phase(9).state(:,2);
states.lat9     = input.phase(9).state(:,3);
states.v9       = input.phase(9).state(:,4);
% states.gamma   = input.phase(9).state(:,5);
% states.zeta    = input.phase(9).state(:,6);
% 
% controls.Alpha = input.phase(9).state(:,7); % Note the 'controls' class here is not the same as the 'control' defined for GPOPS (these are the vehicle controls, not the control theory controls)
% controls.eta   = input.phase(9).state(:,8);
% 
% controls.throttle  = 1;
% 
% Alphadot1  = input.phase(9).control(:,1);
% etadot1 = input.phase(9).control(:,2);
% 
% time = input.phase(9).time;
% 
% auxdata = input.auxdata;

Alphadot9  = input.phase(9).control(:,1)';
etadot9  = input.phase(9).control(:,2)';
time9 = input.phase(9).time;
throttle9 = 0;

[altdot9,londot9,latdot9,gammadot9,vdot9,azidot9, q9, M9, Fd9, rho9,L9,Fueldt9,T9,Isp9,Isp9,m9,heating_rate9,total_acceleration9] = SpaceLinerVehicleModel(time9,input.phase(9),throttle9,auxdata,3);

phaseout(9).dynamics  = [altdot9', londot9', latdot9', vdot9', gammadot9', azidot9',Fueldt9', Alphadot9',etadot9'];


total_acceleration9 = sqrt(vdot9.^2 + (states.v9'.*gammadot9).^2 + (states.v9'.*azidot9).^2)/9.81;

    
phaseout(9).path = [q9,total_acceleration9'];

% phaseout(9).path = [q9, heating_rate9, total_acceleration9']; 

states.lon9(isnan(states.lon9)) = 0;
states.lat9(isnan(states.lat9)) = 0;

states.lon9 = states.lon9 + auxdata.lon0;
states.lon9(states.lon9 > pi) = states.lon9(states.lon9 > pi) - 2*pi;
states.lon9(states.lon9 < -pi) = states.lon9(states.lon9 < -pi) + 2*pi;

% isLand = auxdata.LandSpline(rad2deg(states.lat),rad2deg(states.lon));
% 
% isLand(rad2deg(states.lat) < -60) = 0; % Set antarctica to not be relevant
% 
% AltOverLandCost = zeros(length(states.alt),1);
% AltOverLandCost(isLand==1) = 80000-states.alt(isLand==1); % Find altitudes of every point over land, maximise
% AltOverLandCost(states.alt>80000) = 0;
% % The 'target' must be greater than the max alt, or go to 0 above target,
% % to prevent it being better to fly high over land, than over water
% 
% phaseout(1).integrand = AltOverLandCost.^2/1e6;
% % phaseout(1).integrand = isLand;


AltCost = (80000-states.alt9)/80000;

pop = auxdata.PopInterp(rad2deg(states.lon9),rad2deg(states.lat9));

% popCost = pop.*AltCost; % for flights which go over large amounts of
% population

popCost = pop;

% phaseout(1).integrand = pop;

% Scaling of heating must change if different population density acfuracy
% is used, because population density is per cell, which will change in magnitude.
% phaseout(9).integrand = -(input.phase(1).state(:,1)+input.phase(2).state(:,1)+input.phase(3).state(:,1)+input.phase(4).state(:,1)+input.phase(5).state(:,1)+input.phase(6).state(:,1)+input.phase(7).state(:,1)+input.phase(8).state(:,1)+input.phase(9).state(:,1));
% phaseout(9).integrand = q1+q2+q3+q4+q5+q6+q7+q8+q9;
phaseout(9).integrand = heating_rate9;
% phaseout(1).integrand = pop + heating_rate/1e5;
% phaseout(9).integrand = popCost + heating_rate9/1e5;
end

%======================================================