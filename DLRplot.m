clear all



Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

% auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data
% auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data
% auxdata.interp.T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 
% auxdata.interp.P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 


SLEG = dlmread('SLEG.txt');
TOSCA = dlmread('TOSCA');

HeatMin = dlmread('HeatMin');
HeatMinLimited = dlmread('HeatMinLimited');


t_SLEG = SLEG(:,1);
t_TOSCA = TOSCA(:,1);
t_HeatMin = HeatMin(:,1);
t_HeatMinLimited = HeatMinLimited(:,1);

%%
alt_SLEG = SLEG(:,12);
alt_TOSCA = TOSCA(:,5);
alt_HeatMin = HeatMin(:,5);
alt_HeatMinLimited = HeatMinLimited(:,5);

figure()
hold on
plot(t_SLEG,alt_SLEG/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot(t_TOSCA,alt_TOSCA, 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,alt_HeatMin/1000, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');
plot(t_HeatMinLimited,alt_HeatMinLimited/1000, 'LineWidth', 1, 'color', 'b', 'LineStyle', '-');

xlabel('time (s)');
ylabel('altitude (km)')
legend('SLEG','TOSCA','Base Limits', 'Tight Limits');
%%

v_SLEG = SLEG(:,2);
v_TOSCA = TOSCA(:,2);
v_HeatMin = HeatMin(:,2);
v_HeatMinLimited = HeatMinLimited(:,2);

figure()
hold on
plot(t_SLEG,v_SLEG/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot(t_TOSCA,v_TOSCA, 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,v_HeatMin/1000, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');
plot(t_HeatMinLimited,v_HeatMinLimited/1000, 'LineWidth', 1, 'color', 'b', 'LineStyle', '-');

xlabel('time (s)');
ylabel('Velocity (m/s)')
legend('SLEG','TOSCA','Base Limits', 'Tight Limits');


%%

heat_HeatMin = HeatMin(:,10);
heat_HeatMinLimited = HeatMinLimited(:,10);



rho_SLEG = ppval(rho_spline,alt_SLEG); % Calculate density using atmospheric data
rho_TOSCA = ppval(rho_spline,alt_TOSCA*1e3); % Calculate density using atmospheric data


%Heating model used in Tosca

R_N = 0.205; %effective nose radius (m) 

C = 20254.4;
rho_r = 1.225;
v_r = 10000;
R_Nr = 1;

heat_SLEG = C*sqrt(rho_SLEG/rho_r*R_Nr/R_N).*(v_SLEG/v_r).^3.05*1e4;
heat_TOSCA = C*sqrt(rho_TOSCA/rho_r*R_Nr/R_N).*(v_TOSCA*1000/v_r).^3.05*1e4;


figure()
hold on
plot(t_SLEG,heat_SLEG/1000000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot(t_TOSCA,heat_TOSCA/1000000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,heat_HeatMin/1000000, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');
plot(t_HeatMinLimited,heat_HeatMinLimited/1000000, 'LineWidth', 1, 'color', 'b', 'LineStyle', '-');

xlabel('time (s)');
ylabel('Heating Rate (MW/m^2)')
legend('SLEG','TOSCA','Base Limits', 'Tight Limits');



%%
accel_TOSCA = TOSCA(:,12);
accel_HeatMin = HeatMin(:,8);
accel_HeatMinLimited = HeatMinLimited(:,8);

figure()
hold on
plot(t_TOSCA,accel_TOSCA, 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,accel_HeatMin, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');
plot(t_HeatMinLimited,accel_HeatMinLimited, 'LineWidth', 1, 'color', 'b', 'LineStyle', '-');

xlabel('time (s)');
ylabel('Acceleration (g)')
legend('TOSCA','Base Limits', 'Tight Limits');

%%
q_TOSCA = TOSCA(:,13);
q_HeatMin = HeatMin(:,9);
q_HeatMinLimited = HeatMinLimited(:,9);

q_SLEG = 0.5 * rho_SLEG .* (v_SLEG .^2); % Calculating Dynamic Pressure


figure()
hold on
plot(t_SLEG,q_SLEG/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot(t_TOSCA,q_TOSCA/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,q_HeatMin/1000, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');
plot(t_HeatMinLimited,q_HeatMinLimited/1000, 'LineWidth', 1, 'color', 'b', 'LineStyle', '-');

xlabel('time (s)');
ylabel('Dynamic Pressure (kPa)')
legend('SLEG','TOSCA','Base Limits', 'Tight Limits');