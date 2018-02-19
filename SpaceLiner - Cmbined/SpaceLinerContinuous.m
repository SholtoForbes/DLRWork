function phaseout = SpaceLinerContinuous(input)

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
time1 = input.phase(1).time;
throttle1 = 1;
[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1, M1, Fd1, rho1,L1,Fueldt1,T1,Isp1,Isp2,m1,heating_rate1,total_acceleration1] = SpaceLinerVehicleModel(time1,input.phase(1),throttle1,auxdata,1);
phaseout(1).dynamics  = [altdot1', londot1', latdot1', vdot1', gammadot1', azidot1', -Fueldt1, Alphadot1'];

% phaseout(1).path = [q1,heating_rate1,total_acceleration1'];
phaseout(1).path = [q1,total_acceleration1'];
% phaseout(1).integrand = heating_rate1;
phaseout(1).integrand = q1/1e5;
%% 2
Alphadot2  = input.phase(2).control(:,1)';
time2 = input.phase(2).time;
throttle2 = 0.891;
[altdot2,londot2,latdot2,gammadot2,vdot2,azidot2, q2, M2, Fd2, rho2,L2,Fueldt2,T2,Isp2,Isp2,m2,heating_rate2,total_acceleration2] = SpaceLinerVehicleModel(time2,input.phase(2),throttle2,auxdata,1);
phaseout(2).dynamics  = [altdot2', londot2', latdot2', vdot2', gammadot2', azidot2', -Fueldt2, Alphadot2'];
% total_acceleration2
% phaseout(2).path = [q2,heating_rate2,total_acceleration2'];
phaseout(2).path = [q2,total_acceleration2'];
% phaseout(2).integrand = heating_rate2;
phaseout(2).integrand = q2/1e5;
%% 3
Alphadot3  = input.phase(3).control(:,1)';
time3 = input.phase(3).time;
throttle3 = 0.812;
[altdot3,londot3,latdot3,gammadot3,vdot3,azidot3, q3, M3, Fd3, rho3,L3,Fueldt3,T3,Isp3,Isp3,m3,heating_rate3,total_acceleration3] = SpaceLinerVehicleModel(time3,input.phase(3),throttle3,auxdata,1);
phaseout(3).dynamics  = [altdot3', londot3', latdot3', vdot3', gammadot3', azidot3', -Fueldt3, Alphadot3'];
% phaseout(3).path = [q3,heating_rate3,total_acceleration3'];
phaseout(3).path = [q3,total_acceleration3'];
% phaseout(3).integrand = heating_rate3;
phaseout(3).integrand = q3/1e4;
%% 4
Alphadot4  = input.phase(4).control(:,1)';
time4 = input.phase(4).time;
throttle4 = .7333;
[altdot4,londot4,latdot4,gammadot4,vdot4,azidot4, q4, M4, Fd4, rho4,L4,Fueldt4,T4,Isp4,Isp4,m4,heating_rate4,total_acceleration4] = SpaceLinerVehicleModel(time4,input.phase(4),throttle4,auxdata,1);
phaseout(4).dynamics  = [altdot4', londot4', latdot4', vdot4', gammadot4', azidot4', -Fueldt4, Alphadot4'];
% phaseout(4).path = [q4,heating_rate4,total_acceleration4'];
phaseout(4).path = [q4,total_acceleration4'];
% phaseout(4).integrand = heating_rate4;
phaseout(4).integrand = q4/1e5;
%% 5
Alphadot5  = input.phase(5).control(:,1)';
time5 = input.phase(5).time;
throttle5 = .6545;
[altdot5,londot5,latdot5,gammadot5,vdot5,azidot5, q5, M5, Fd5, rho5,L5,Fueldt5,T5,Isp5,Isp5,m5,heating_rate5,total_acceleration5] = SpaceLinerVehicleModel(time5,input.phase(5),throttle5,auxdata,1);
phaseout(5).dynamics  = [altdot5', londot5', latdot5', vdot5', gammadot5', azidot5', -Fueldt5, Alphadot5'];
% phaseout(5).path = [q5,heating_rate5,total_acceleration5'];
phaseout(5).path = [q5,total_acceleration5'];
% phaseout(5).integrand = heating_rate5;
phaseout(5).integrand = q5/1e5;
%% 6
Alphadot6  = input.phase(6).control(:,1)';
time6 = input.phase(6).time;
throttle6 = .5757;
[altdot6,londot6,latdot6,gammadot6,vdot6,azidot6, q6, M6, Fd6, rho6,L6,Fueldt6,T6,Isp6,Isp6,m6,heating_rate6,total_acceleration6] = SpaceLinerVehicleModel(time6,input.phase(6),throttle6,auxdata,1);
phaseout(6).dynamics  = [altdot6', londot6', latdot6', vdot6', gammadot6', azidot6', -Fueldt6, Alphadot6'];
% phaseout(6).path = [q6,heating_rate6,total_acceleration6'];
phaseout(6).path = [q6,total_acceleration6'];
% phaseout(6).integrand = heating_rate6;
phaseout(6).integrand = q6/1e5;
%% 7
Alphadot7  = input.phase(7).control(:,1)';
time7 = input.phase(7).time;
throttle7 = 0.496;
[altdot7,londot7,latdot7,gammadot7,vdot7,azidot7, q7, M7, Fd7, rho7,L7,Fueldt7,T7,Isp7,Isp7,m7,heating_rate7,total_acceleration7] = SpaceLinerVehicleModel(time7,input.phase(7),throttle7,auxdata,1);
phaseout(7).dynamics  = [altdot7', londot7', latdot7', vdot7', gammadot7', azidot7', -Fueldt7, Alphadot7'];
% phaseout(7).path = [q7,heating_rate7,total_acceleration7'];
phaseout(7).path = [q7,total_acceleration7'];
% phaseout(7).integrand = heating_rate7;
phaseout(7).integrand = q7/1e5;
%% 8
Alphadot8  = input.phase(8).control(:,1)';
time8 = input.phase(8).time;
throttle8 = 1;
[altdot8,londot8,latdot8,gammadot8,vdot8,azidot8, q8, M8, Fd8, rho8,L8,Fueldt8,T8,Isp8,Isp8,m8,heating_rate8,total_acceleration8] = SpaceLinerVehicleModel(time8,input.phase(8),throttle8,auxdata,2);
phaseout(8).dynamics  = [altdot8', londot8', latdot8', vdot8', gammadot8', azidot8', -Fueldt8, Alphadot8'];
% phaseout(8).path = [q8,heating_rate8,total_acceleration8'];
phaseout(8).path = [q8,total_acceleration8'];
% phaseout(8).integrand = heating_rate8;
phaseout(8).integrand = q8/1e5;
%% 9
Alphadot9  = input.phase(9).control(:,1)';
time9 = input.phase(9).time;
throttle9 = 0;
[altdot9,londot9,latdot9,gammadot9,vdot9,azidot9, q9, M9, Fd9, rho9,L9,Fueldt9,T9,Isp9,Isp9,m9,heating_rate9,total_acceleration9] = SpaceLinerVehicleModel(time9,input.phase(9),throttle9,auxdata,2);
phaseout(9).dynamics  = [altdot9', londot9', latdot9', vdot9', gammadot9', azidot9', -Fueldt9, Alphadot9'];
% phaseout(9).path = [q9,heating_rate9,total_acceleration9'];
phaseout(9).path = [q9,total_acceleration9'];
% phaseout(9).integrand = heating_rate9;
phaseout(9).integrand = q9/1e5;
end

%======================================================