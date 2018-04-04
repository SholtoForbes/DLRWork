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
etadot1  = input.phase(1).control(:,2)';
time1 = input.phase(1).time;
throttle1 = 1;
[altdot1,londot1,latdot1,gammadot1,vdot1,azidot1, q1, M1, Fd1, rho1,L1,Fueldt1,T1,Isp1,Isp2,m1,heating_rate1,total_acceleration1] = SpaceLinerVehicleModel(time1,input.phase(1),throttle1,auxdata,1);
phaseout(1).dynamics  = [altdot1', londot1', latdot1', vdot1', gammadot1', azidot1', -Fueldt1, Alphadot1',etadot1'];

% phaseout(1).path = [q1,total_acceleration1'];
phaseout(1).path = [q1, heating_rate1,total_acceleration1'];
% phaseout(1).integrand = heating_rate1;

%
alt1     = input.phase(1).state(:,1);
lon1     = input.phase(1).state(:,2);
lat1     = input.phase(1).state(:,3);

% lon1(isnan(lon1)) = 0;
% lat1(isnan(lat1)) = 0;

popCost1 = ones(length(lon1),1);
popCost1(isnan(lon1)) = nan;
popCost1(isnan(lat1)) = nan;
lon1 = lon1 + auxdata.lon0;
lon1(lon1 > pi) = lon1(lon1 > pi) - 2*pi;
lon1(lon1 < -pi) = lon1(lon1 < -pi) + 2*pi;

AltCost1 = (80000-alt1)/10000;
AltCost1(alt1>80000) = 0;
pop1 = NaN*ones(length(lon1),1);
pop1(not(isnan(popCost1))) = auxdata.PopInterp(rad2deg(lon1(not(isnan(popCost1)))),rad2deg(lat1(not(isnan(popCost1)))));
popCost1 = pop1.*AltCost1; % for flights which go over large amounts of
% popCost1 = popCost1.*L1/10^4;
% heating_rate_cost1 = heating_rate1;
% heating_rate_cost1(heating_rate1<1.3e6) = 0;
% phaseout(1).integrand =  heating_rate_cost1;

phaseout(1).integrand = popCost1 + heating_rate1/1e5;
% phaseout(1).integrand = -alt1;
% phaseout(1).integrand = heating_rate1/1e5;
%% 2
Alphadot2  = input.phase(2).control(:,1)';
etadot2  = input.phase(2).control(:,2)';
time2 = input.phase(2).time;
throttle2 = 0.891;
[altdot2,londot2,latdot2,gammadot2,vdot2,azidot2, q2, M2, Fd2, rho2,L2,Fueldt2,T2,Isp2,Isp2,m2,heating_rate2,total_acceleration2] = SpaceLinerVehicleModel(time2,input.phase(2),throttle2,auxdata,1);
phaseout(2).dynamics  = [altdot2', londot2', latdot2', vdot2', gammadot2', azidot2', -Fueldt2, Alphadot2',etadot2'];
% total_acceleration2
% phaseout(2).path = [q2,total_acceleration2'];
phaseout(2).path = [q2,heating_rate2,total_acceleration2'];
% phaseout(2).integrand = heating_rate2;

%
alt2     = input.phase(2).state(:,1);
lon2     = input.phase(2).state(:,2);
lat2     = input.phase(2).state(:,3);

% lon2(isnan(lon2)) = 0;
% lat2(isnan(lat2)) = 0;

% if any(isnan(lon2)) || any(isnan(lat2))

popCost2 = ones(length(lon2),1);
popCost2(isnan(lon2)) = nan;
popCost2(isnan(lat2)) = nan;

lon2 = lon2 + auxdata.lon0;

lon2(lon2 > pi) = lon2(lon2 > pi) - 2*pi;
lon2(lon2 < -pi) = lon2(lon2 < -pi) + 2*pi;

AltCost2 = (80000-alt2)/10000;
AltCost2(alt2>80000) = 0;
pop2 = NaN*ones(length(lon2),1);
pop2(not(isnan(popCost2))) = auxdata.PopInterp(rad2deg(lon2(not(isnan(popCost2)))),rad2deg(lat2(not(isnan(popCost2)))));
popCost2 = pop2.*AltCost2; % for flights which go over large amounts of

% popCost2 = popCost2.*L2/10^4;

% end
% phaseout(2).path = [q2,heating_rate2,total_acceleration2',AltCost2];

% phaseout(2).integrand = heating_rate2/1e5;
% heating_rate_cost2(heating_rate2<1.3e6) = 0;
% phaseout(2).integrand = popCost2 + heating_rate_cost2;

phaseout(2).integrand = popCost2 + heating_rate2/1e5;
% phaseout(2).integrand = -alt2;
%% 3
Alphadot3  = input.phase(3).control(:,1)';
etadot3  = input.phase(3).control(:,2)';
time3 = input.phase(3).time;
throttle3 = 0.812;
[altdot3,londot3,latdot3,gammadot3,vdot3,azidot3, q3, M3, Fd3, rho3,L3,Fueldt3,T3,Isp3,Isp3,m3,heating_rate3,total_acceleration3] = SpaceLinerVehicleModel(time3,input.phase(3),throttle3,auxdata,1);
phaseout(3).dynamics  = [altdot3', londot3', latdot3', vdot3', gammadot3', azidot3', -Fueldt3, Alphadot3',etadot3'];
% phaseout(3).path = [q3,total_acceleration3'];
phaseout(3).path = [q3,heating_rate3,total_acceleration3'];
% phaseout(3).integrand = heating_rate3;

%
alt3     = input.phase(3).state(:,1);
lon3     = input.phase(3).state(:,2);
lat3     = input.phase(3).state(:,3);

% lon3(isnan(lon3)) = 0;
% lat3(isnan(lat3)) = 0;

popCost3 = ones(length(lon3),1);
popCost3(isnan(lon3)) = nan;
popCost3(isnan(lat3)) = nan;
lon3 = lon3 + auxdata.lon0;
lon3(lon3 > pi) = lon3(lon3 > pi) - 2*pi;
lon3(lon3 < -pi) = lon3(lon3 < -pi) + 2*pi;

AltCost3 = (80000-alt3)/10000;
AltCost3(alt3>80000) = 0;
pop3 = NaN*ones(length(lon3),1);
pop3(not(isnan(popCost3))) = auxdata.PopInterp(rad2deg(lon3(not(isnan(popCost3)))),rad2deg(lat3(not(isnan(popCost3)))));
popCost3 = pop3.*AltCost3; % for flights which go over large amounts of

% popCost3 = popCost3.*L3/10^4;
% phaseout(3).path = [q3,heating_rate3,total_acceleration3',AltCost3];

% phaseout(3).integrand = heating_rate3/1e5;
% heating_rate_cost3(heating_rate3<1.3e6) = 0;
% phaseout(3).integrand = popCost3 + heating_rate_cost3;
phaseout(3).integrand = popCost3 + heating_rate3/1e5;
% phaseout(3).integrand = -alt3;
%% 4
Alphadot4  = input.phase(4).control(:,1)';
etadot4  = input.phase(4).control(:,2)';
time4 = input.phase(4).time;
throttle4 = .7333;
[altdot4,londot4,latdot4,gammadot4,vdot4,azidot4, q4, M4, Fd4, rho4,L4,Fueldt4,T4,Isp4,Isp4,m4,heating_rate4,total_acceleration4] = SpaceLinerVehicleModel(time4,input.phase(4),throttle4,auxdata,1);
phaseout(4).dynamics  = [altdot4', londot4', latdot4', vdot4', gammadot4', azidot4', -Fueldt4, Alphadot4',etadot4'];
% phaseout(4).path = [q4,total_acceleration4'];
phaseout(4).path = [q4,heating_rate4,total_acceleration4'];
% phaseout(4).integrand = heating_rate4;

%
alt4     = input.phase(4).state(:,1);
lon4     = input.phase(4).state(:,2);
lat4     = input.phase(4).state(:,3);

% lon4(isnan(lon4)) = 0;
% lat4(isnan(lat4)) = 0;
popCost4 = ones(length(lon4),1);
popCost4(isnan(lon4)) = nan;
popCost4(isnan(lat4)) = nan;
lon4 = lon4 + auxdata.lon0;
lon4(lon4 > pi) = lon4(lon4 > pi) - 2*pi;
lon4(lon4 < -pi) = lon4(lon4 < -pi) + 2*pi;

AltCost4 = (80000-alt4)/10000;
AltCost4(alt4>80000) = 0;
pop4 = NaN*ones(length(lon4),1);
pop4(not(isnan(popCost4))) = auxdata.PopInterp(rad2deg(lon4(not(isnan(popCost4)))),rad2deg(lat4(not(isnan(popCost4)))));
popCost4 = pop4.*AltCost4; % for flights which go over large amounts of
% popCost4 = popCost4.*L4/11^4;
% phaseout(4).path = [q4,heating_rate4,total_acceleration4',AltCost4];

% phaseout(4).integrand = heating_rate4/1e5;
% heating_rate_cost4(heating_rate4<1.3e6) = 0;
% phaseout(4).integrand = popCost4 + heating_rate_cost4;
phaseout(4).integrand = popCost4 + heating_rate4/1e5;
% phaseout(4).integrand = -alt4;
%% 5
Alphadot5  = input.phase(5).control(:,1)';
etadot5  = input.phase(5).control(:,2)';
time5 = input.phase(5).time;
throttle5 = .6545;
[altdot5,londot5,latdot5,gammadot5,vdot5,azidot5, q5, M5, Fd5, rho5,L5,Fueldt5,T5,Isp5,Isp5,m5,heating_rate5,total_acceleration5] = SpaceLinerVehicleModel(time5,input.phase(5),throttle5,auxdata,1);
phaseout(5).dynamics  = [altdot5', londot5', latdot5', vdot5', gammadot5', azidot5', -Fueldt5, Alphadot5',etadot5'];
% phaseout(5).path = [q5,total_acceleration5'];
phaseout(5).path = [q5,heating_rate5,total_acceleration5'];
% phaseout(5).integrand = heating_rate5;

%
alt5     = input.phase(5).state(:,1);
lon5     = input.phase(5).state(:,2);
lat5     = input.phase(5).state(:,3);

% lon5(isnan(lon5)) = 0;
% lat5(isnan(lat5)) = 0;
popCost5 = ones(length(lon5),1);
popCost5(isnan(lon5)) = nan;
popCost5(isnan(lat5)) = nan;
lon5 = lon5 + auxdata.lon0;
lon5(lon5 > pi) = lon5(lon5 > pi) - 2*pi;
lon5(lon5 < -pi) = lon5(lon5 < -pi) + 2*pi;

AltCost5 = (80000-alt5)/10000;
AltCost5(alt5>80000) = 0;
pop5 = NaN*ones(length(lon5),1);
pop5(not(isnan(popCost5))) = auxdata.PopInterp(rad2deg(lon5(not(isnan(popCost5)))),rad2deg(lat5(not(isnan(popCost5)))));
popCost5 = pop5.*AltCost5; % for flights which go over large amounts of
% popCost5 = popCost5.*L5/10^4;
% phaseout(5).path = [q5,heating_rate5,total_acceleration5',AltCost5];

% phaseout(5).integrand = heating_rate5/1e6;
% heating_rate_cost5(heating_rate5<1.3e6) = 0;
% phaseout(5).integrand = popCost5 + heating_rate_cost5;

phaseout(5).integrand = popCost5 + heating_rate5/1e5;
% phaseout(5).integrand = -alt5;
%% 6
Alphadot6  = input.phase(6).control(:,1)';
etadot6  = input.phase(6).control(:,2)';
time6 = input.phase(6).time;
throttle6 = .5757;
[altdot6,londot6,latdot6,gammadot6,vdot6,azidot6, q6, M6, Fd6, rho6,L6,Fueldt6,T6,Isp6,Isp6,m6,heating_rate6,total_acceleration6] = SpaceLinerVehicleModel(time6,input.phase(6),throttle6,auxdata,1);
phaseout(6).dynamics  = [altdot6', londot6', latdot6', vdot6', gammadot6', azidot6', -Fueldt6, Alphadot6',etadot6'];
% phaseout(6).path = [q6,total_acceleration6'];
phaseout(6).path = [q6,heating_rate6,total_acceleration6'];
% phaseout(6).integrand = heating_rate6;

%
alt6     = input.phase(6).state(:,1);
lon6     = input.phase(6).state(:,2);
lat6     = input.phase(6).state(:,3);

% lon6(isnan(lon6)) = 0;
% lat6(isnan(lat6)) = 0;
popCost6 = ones(length(lon6),1);
popCost6(isnan(lon6)) = nan;
popCost6(isnan(lat6)) = nan;
lon6 = lon6 + auxdata.lon0;
lon6(lon6 > pi) = lon6(lon6 > pi) - 2*pi;
lon6(lon6 < -pi) = lon6(lon6 < -pi) + 2*pi;

AltCost6 = (80000-alt6)/10000;
AltCost6(alt6>80000) = 0;
pop6 = NaN*ones(length(lon6),1);
pop6(not(isnan(popCost6))) = auxdata.PopInterp(rad2deg(lon6(not(isnan(popCost6)))),rad2deg(lat6(not(isnan(popCost6)))));
popCost6 = pop6.*AltCost6; % for flights which go over large amounts of
% popCost6 = popCost6.*L6/10^4;
% phaseout(6).path = [q6,heating_rate6,total_acceleration6',AltCost6];
% phaseout(6).integrand = heating_rate6/1e5;
% heating_rate_cost6(heating_rate6<1.3e6) = 0;
% phaseout(6).integrand = popCost6 + heating_rate_cost6;
phaseout(6).integrand = popCost6 + heating_rate6/1e5;
% phaseout(6).integrand = -alt6;
%% 7
Alphadot7  = input.phase(7).control(:,1)';
etadot7  = input.phase(7).control(:,2)';
time7 = input.phase(7).time;
throttle7 = 0.496;
[altdot7,londot7,latdot7,gammadot7,vdot7,azidot7, q7, M7, Fd7, rho7,L7,Fueldt7,T7,Isp7,Isp7,m7,heating_rate7,total_acceleration7] = SpaceLinerVehicleModel(time7,input.phase(7),throttle7,auxdata,1);
phaseout(7).dynamics  = [altdot7', londot7', latdot7', vdot7', gammadot7', azidot7', -Fueldt7, Alphadot7',etadot7'];
% phaseout(7).path = [q7,total_acceleration7'];
phaseout(7).path = [q7,heating_rate7,total_acceleration7'];
% phaseout(7).integrand = heating_rate7;

%
alt7     = input.phase(7).state(:,1);
lon7     = input.phase(7).state(:,2);
lat7     = input.phase(7).state(:,3);

% lon7(isnan(lon7)) = 0;
% lat7(isnan(lat7)) = 0;
popCost7 = ones(length(lon7),1);
popCost7(isnan(lon7)) = nan;
popCost7(isnan(lat7)) = nan;
lon7 = lon7 + auxdata.lon0;
lon7(lon7 > pi) = lon7(lon7 > pi) - 2*pi;
lon7(lon7 < -pi) = lon7(lon7 < -pi) + 2*pi;

AltCost7 = (80000-alt7)/10000;
AltCost7(alt7>80000) = 0;
pop7 = NaN*ones(length(lon7),1);
pop7(not(isnan(popCost7))) = auxdata.PopInterp(rad2deg(lon7(not(isnan(popCost7)))),rad2deg(lat7(not(isnan(popCost7)))));
popCost7 = pop7.*AltCost7; % for flights which go over large amounts of
% popCost7 = popCost7.*L7/10^5;
% phaseout(7).path = [q7,heating_rate7,total_acceleration7',AltCost7];
% phaseout(7).integrand = heating_rate7/1e5;
% heating_rate_cost7(heating_rate7<1.3e6) = 0;
% phaseout(7).integrand = popCost7 + heating_rate_cost7;
phaseout(7).integrand = popCost7 + heating_rate7/1e5;
% phaseout(7).integrand = -alt7;
%% 8
Alphadot8  = input.phase(8).control(:,1)';
etadot8  = input.phase(8).control(:,2)';
time8 = input.phase(8).time;
throttle8 = 1;
[altdot8,londot8,latdot8,gammadot8,vdot8,azidot8, q8, M8, Fd8, rho8,L8,Fueldt8,T8,Isp8,Isp8,m8,heating_rate8,total_acceleration8] = SpaceLinerVehicleModel(time8,input.phase(8),throttle8,auxdata,2);
phaseout(8).dynamics  = [altdot8', londot8', latdot8', vdot8', gammadot8', azidot8', -Fueldt8, Alphadot8',etadot8'];
% phaseout(8).path = [q8,total_acceleration8'];
phaseout(8).path = [q8,heating_rate8,total_acceleration8'];
% phaseout(8).integrand = heating_rate8;

%
alt8     = input.phase(8).state(:,1);
lon8     = input.phase(8).state(:,2);
lat8     = input.phase(8).state(:,3);

% lon8(isnan(lon8)) = 0;
% lat8(isnan(lat8)) = 0;
popCost8 = ones(length(lon8),1);
popCost8(isnan(lon8)) = nan;
popCost8(isnan(lat8)) = nan;
lon8 = lon8 + auxdata.lon0;
lon8(lon8 > pi) = lon8(lon8 > pi) - 2*pi;
lon8(lon8 < -pi) = lon8(lon8 < -pi) + 2*pi;

AltCost8 = (80000-alt8)/10000;
AltCost8(alt8>80000) = 0;
pop8 = NaN*ones(length(lon8),1);
pop8(not(isnan(popCost8))) = auxdata.PopInterp(rad2deg(lon8(not(isnan(popCost8)))),rad2deg(lat8(not(isnan(popCost8)))));
popCost8 = pop8.*AltCost8; % for flights which go over large amounts of
% popCost8 = popCost8.*L8/10^4;
% phaseout(8).path = [q8,heating_rate8,total_acceleration8',AltCost8];
% phaseout(8).integrand = heating_rate8/1e5;
% heating_rate_cost8(heating_rate8<1.3e6) = 0;
% phaseout(8).integrand = popCost8 + heating_rate_cost8;
phaseout(8).integrand = popCost8 + heating_rate8/1e5;
% phaseout(8).integrand = -alt8;
%% 9
Alphadot9  = input.phase(9).control(:,1)';
etadot9  = input.phase(9).control(:,2)';
time9 = input.phase(9).time;
throttle9 = 0.86667;
[altdot9,londot9,latdot9,gammadot9,vdot9,azidot9, q9, M9, Fd9, rho9,L9,Fueldt9,T9,Isp9,Isp9,m9,heating_rate9,total_acceleration9] = SpaceLinerVehicleModel(time9,input.phase(9),throttle9,auxdata,2);
phaseout(9).dynamics  = [altdot9', londot9', latdot9', vdot9', gammadot9', azidot9', -Fueldt9, Alphadot9',etadot9'];
% phaseout(9).path = [q9,total_acceleration9'];
phaseout(9).path = [q9,heating_rate9,total_acceleration9'];
% phaseout(9).integrand = heating_rate9;

%
alt9     = input.phase(9).state(:,1);
lon9     = input.phase(9).state(:,2);
lat9     = input.phase(9).state(:,3);

% lon9(isnan(lon9)) = 0;
% lat9(isnan(lat9)) = 0;
popCost9 = ones(length(lon9),1);
popCost9(isnan(lon9)) = nan;
popCost9(isnan(lat9)) = nan;
lon9 = lon9 + auxdata.lon0;
lon9(lon9 > pi) = lon9(lon9 > pi) - 2*pi;
lon9(lon9 < -pi) = lon9(lon9 < -pi) + 2*pi;

AltCost9 = (80000-alt9)/10000;
AltCost9(alt9>80000) = 0;
pop9 = NaN*ones(length(lon9),1);
pop9(not(isnan(popCost9))) = auxdata.PopInterp(rad2deg(lon9(not(isnan(popCost9)))),rad2deg(lat9(not(isnan(popCost9)))));
popCost9 = pop9.*AltCost9; % for flights which go over large amounts of
% popCost9 = popCost9.*L9/10^4;
% phaseout(9).path = [q9,heating_rate9,total_acceleration9',AltCost9];
% phaseout(9).integrand = heating_rate9/1e5;
% heating_rate_cost9(heating_rate9<1.3e6) = 0;
% phaseout(9).integrand = popCost9 + heating_rate_cost9;
phaseout(9).integrand = popCost9 + heating_rate9/1e5;
% phaseout(9).integrand = -alt9;
%% 10
Alphadot10  = input.phase(10).control(:,1)';
etadot10  = input.phase(10).control(:,2)';
time10 = input.phase(10).time;
throttle10 = 0;
[altdot10,londot10,latdot10,gammadot10,vdot10,azidot10, q10, M10, Fd10, rho10,L10,Fueldt10,T10,Isp10,Isp10,m10,heating_rate10,total_acceleration10] = SpaceLinerVehicleModel(time10,input.phase(10),throttle10,auxdata,2);
phaseout(10).dynamics  = [altdot10', londot10', latdot10', vdot10', gammadot10', azidot10', -Fueldt10, Alphadot10',etadot10'];
% phaseout(10).path = [q10,total_acceleration10'];
phaseout(10).path = [q10,heating_rate10,total_acceleration10'];
% phaseout(10).integrand = heating_rate10;

%
alt10     = input.phase(10).state(:,1);
lon10     = input.phase(10).state(:,2);
lat10     = input.phase(10).state(:,3);

% lon10(isnan(lon10)) = 0;
% lat10(isnan(lat10)) = 0;
popCost10 = ones(length(lon10),1);
popCost10(isnan(lon10)) = nan;
popCost10(isnan(lat10)) = nan;
lon10 = lon10 + auxdata.lon0;
lon10(lon10 > pi) = lon10(lon10 > pi) - 2*pi;
lon10(lon10 < -pi) = lon10(lon10 < -pi) + 2*pi;

AltCost10 = (80000-alt10)/10000;
AltCost10(alt10>80000) = 0;
pop10 = NaN*ones(length(lon10),1);
pop10(not(isnan(popCost10))) = auxdata.PopInterp(rad2deg(lon10(not(isnan(popCost10)))),rad2deg(lat10(not(isnan(popCost10)))));
popCost10 = pop10.*AltCost10; % for flights which go over large amounts of
% popCost10 = popCost10.*L10/10^4;
% phaseout(10).integrand = heating_rate10/1e5;
% heating_rate_cost10(heating_rate10<1.3e6) = 0;
% phaseout(10).integrand = popCost10 + heating_rate_cost10;
phaseout(10).integrand = popCost10 + heating_rate10/1e5;
% phaseout(10).integrand = -alt10;
end

%======================================================