clear all



Atmosphere = dlmread('atmosphere.txt');
interp.Atmosphere = Atmosphere;
auxdata.interp.Atmosphere = interp.Atmosphere;

% auxdata.interp.c_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,5)); % Calculate speed of sound using atmospheric data
rho_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,4)); % Calculate density using atmospheric data
% auxdata.interp.p_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); % Calculate density using atmospheric data
% auxdata.interp.T0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,2)); 
% auxdata.interp.P0_spline = spline( interp.Atmosphere(:,1),  interp.Atmosphere(:,3)); 

TOSCA_ascent = dlmread('tosca_ascent.data');

SLEG = dlmread('SLEG.txt');
TOSCA = dlmread('TOSCA');

out = dlmread('out');

t_TOSCA_ascent = TOSCA_ascent(:,1);

t_SLEG = SLEG(:,1);
t_TOSCA = TOSCA(:,1);
t_HeatMin = out(:,1);


%%

alt_TOSCA_ascent = TOSCA_ascent(:,5);

alt_SLEG = SLEG(:,12);
alt_TOSCA = TOSCA(:,5);
alt_HeatMin = out(:,5);

figure()
hold on
% plot(t_TOSCA_ascent,alt_TOSCA_ascent, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');

plot(t_SLEG+t_TOSCA_ascent(end),alt_SLEG/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot([t_TOSCA_ascent' t_TOSCA'+t_TOSCA_ascent(end)],[alt_TOSCA_ascent' alt_TOSCA'], 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,alt_HeatMin, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');

xlabel('time (s)');
ylabel('altitude (km)')
legend('SLEG','TOSCA', 'Optimised');
%%
v_TOSCA_ascent = TOSCA_ascent(:,2);

v_SLEG = SLEG(:,2);
v_TOSCA = TOSCA(:,2);
v_HeatMin = out(:,2);

figure()
hold on
% plot(t_TOSCA_ascent,v_TOSCA_ascent, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');

plot(t_SLEG+t_TOSCA_ascent(end),v_SLEG/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot([t_TOSCA_ascent' t_TOSCA'+t_TOSCA_ascent(end)],[v_TOSCA_ascent' v_TOSCA'], 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,v_HeatMin, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');

xlabel('time (s)');
ylabel('Velocity (m/s)')
legend('SLEG','TOSCA', 'Optimised');


%%



heat_HeatMin = out(:,14);

rho_TOSCA_ascent = ppval(rho_spline,alt_TOSCA_ascent*1e3); % Calculate density using atmospheric data

rho_SLEG = ppval(rho_spline,alt_SLEG); % Calculate density using atmospheric data
rho_TOSCA = ppval(rho_spline,alt_TOSCA*1e3); % Calculate density using atmospheric data


%Heating model used in Tosca

R_N = 0.205; %effective nose radius (m) 

C = 20254.4;
rho_r = 1.225;
v_r = 10000;
R_Nr = 1;
heat_TOSCA_ascent = C*sqrt(rho_TOSCA_ascent/rho_r*R_Nr/R_N).*(v_TOSCA_ascent*1000/v_r).^3.05*1e4;

heat_SLEG = C*sqrt(rho_SLEG/rho_r*R_Nr/R_N).*(v_SLEG/v_r).^3.05*1e4;
heat_TOSCA = C*sqrt(rho_TOSCA/rho_r*R_Nr/R_N).*(v_TOSCA*1000/v_r).^3.05*1e4;


figure()
hold on
% plot(t_TOSCA_ascent,heat_TOSCA_ascent/1000000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');

plot(t_SLEG+t_TOSCA_ascent(end),heat_SLEG/1000000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot([t_TOSCA_ascent' t_TOSCA'+t_TOSCA_ascent(end)],[heat_TOSCA_ascent' heat_TOSCA']/1000000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,heat_HeatMin/1000000, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');

xlabel('time (s)');
ylabel('Heating Rate (MW/m^2)')
legend('SLEG','TOSCA', 'Optimised');



%%
accel_TOSCA_ascent = TOSCA_ascent(:,12);

accel_TOSCA = TOSCA(:,12);
accel_HeatMin = out(:,12);

figure()
hold on
% plot(t_TOSCA_ascent,accel_TOSCA_ascent, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');

plot([t_TOSCA_ascent' t_TOSCA'+t_TOSCA_ascent(end)],[accel_TOSCA_ascent' accel_TOSCA'], 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,accel_HeatMin/9.81, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');

xlabel('time (s)');
ylabel('Acceleration (g)')
legend('TOSCA', 'Optimised');

%%
q_TOSCA_ascent = TOSCA_ascent(:,13);
q_TOSCA = TOSCA(:,13);
q_HeatMin = out(:,13);

q_SLEG = 0.5 * rho_SLEG .* (v_SLEG .^2); % Calculating Dynamic Pressure


figure()
hold on
% plot(t_TOSCA_ascent,q_TOSCA_ascent/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');

plot(t_SLEG+t_TOSCA_ascent(end),q_SLEG/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');
plot([t_TOSCA_ascent' t_TOSCA'+t_TOSCA_ascent(end)],[q_TOSCA_ascent' q_TOSCA']/1000, 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(t_HeatMin,q_HeatMin/1000, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');

xlabel('time (s)');
ylabel('Dynamic Pressure (kPa)')
legend('SLEG','TOSCA', 'Optimised');


%%

lat_TOSCA_ascent = TOSCA_ascent(:,7);
lat_TOSCA = TOSCA(:,7);
lat_HeatMin = out(:,7);

lon_TOSCA_ascent = TOSCA_ascent(:,6);
lon_TOSCA = TOSCA(:,6);
lon_HeatMin = out(:,6);

figure()
hold on
% plot(t_TOSCA_ascent,alt_TOSCA_ascent, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');

plot([lon_TOSCA_ascent' lon_TOSCA'],[lat_TOSCA_ascent' lat_TOSCA'], 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(lon_HeatMin,lat_HeatMin, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');

xlabel('Longitude');
ylabel('Latitude')
legend('TOSCA', 'Optimised');