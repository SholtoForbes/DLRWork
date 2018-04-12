clear all

lon0 = 150.5100; % Rockhampton 

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

Q_SLEG = cumtrapz([t_TOSCA_ascent' t_SLEG'+t_TOSCA_ascent(end)],[heat_TOSCA_ascent' heat_SLEG']);
Q_TOSCA = cumtrapz([t_TOSCA_ascent' t_TOSCA'+t_TOSCA_ascent(end)],[heat_TOSCA_ascent' heat_TOSCA']);
Q_HeatMin = cumtrapz(t_HeatMin,heat_HeatMin);

figure()
hold on
plot([t_TOSCA_ascent' t_SLEG'+t_TOSCA_ascent(end)],Q_SLEG);
plot(t_HeatMin,Q_HeatMin);


disp(['SLEG Integrated Heat Load Decreased By ' num2str((1-Q_HeatMin(end)/Q_SLEG(end))*100) '%'])
disp(['TOSCA Integrated Heat Load Decreased By ' num2str((1-Q_HeatMin(end)/Q_TOSCA(end))*100) '%'])

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

lat_SLEG = rad2deg(SLEG(:,7));

lon_TOSCA_ascent = TOSCA_ascent(:,6);
lon_TOSCA = TOSCA(:,6);
lon_HeatMin = out(:,6);


lon_SLEG = rad2deg(SLEG(:,6));


figure()
hold on
% plot(t_TOSCA_ascent,alt_TOSCA_ascent, 'LineWidth', 1, 'color', 'k', 'LineStyle', '-');

plot([lon_TOSCA_ascent' lon_TOSCA'],[lat_TOSCA_ascent' lat_TOSCA'], 'LineWidth', 1, 'color', 'k', 'LineStyle', '--');
plot(lon_HeatMin+lon0,lat_HeatMin, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--');

xlabel('Longitude');
ylabel('Latitude')
legend('TOSCA', 'Optimised');



%%

    figure(231)
hold on

axesm('pcarree','Origin',[0 rad2deg(lon0)-100 0])
geoshow('landareas.shp','FaceColor',[0.8 .8 0.8])
plotm([lat_TOSCA_ascent' lat_TOSCA'],[lon_TOSCA_ascent' lon_TOSCA'], 'LineWidth', 2, 'color', 'k', 'LineStyle', '--')
plotm([lat_TOSCA_ascent' lat_SLEG'],[lon_TOSCA_ascent' lon_SLEG'], 'LineWidth', 2)
plotm(lat_HeatMin,lon_HeatMin, 'LineWidth', 0.8, 'color', 'r', 'LineStyle', '--', 'LineWidth', 2)
    cities = shaperead('worldcities', 'UseGeoCoords', true);
lats = extractfield(cities,'Lat');
lons = extractfield(cities,'Lon');
geoshow(lats, lons,...
        'DisplayType', 'point',...
        'Marker', 'o',...
        'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r',...
        'MarkerSize', 2)
legend('TOSCA','SLEG', 'Optimised');
% 



%%


load PopInterp

auxdata.PopInterp = PopInterp;

%
altpop     = alt_HeatMin*1000;
lonpop     = lon_HeatMin;
latpop     = lat_HeatMin;

% lonpop = lonpop + lon0;

lonpop(lonpop > 180) = lonpop(lonpop > 180) - 360;
lonpop(lonpop < -180) = lonpop(lonpop < -180) + 360;

AltCost1 = (80000-altpop)/10000;
AltCost1(alt_HeatMin*1000>80000) = 0;
% AltCost1(altpop>90000) = gaussmf(altpop(altpop>90000),[10000 90000]);

pop1 = auxdata.PopInterp(lonpop,latpop);

popCost1 = pop1.*AltCost1; % for flights which go over large amounts of



altpop     = [alt_TOSCA_ascent'*1000 alt_SLEG'];
lonpop     = [lon_TOSCA_ascent' lon_SLEG'];
latpop     = [lat_TOSCA_ascent' lat_SLEG'];

lonpop(lonpop > 180) = lonpop(lonpop > 180) - 360;
lonpop(lonpop < -180) = lonpop(lonpop < -180) + 360;

AltCost2 = (80000-altpop)/10000;
AltCost2(alt_HeatMin*1000>80000) = 0;
% AltCost1(altpop>90000) = gaussmf(altpop(altpop>90000),[10000 90000]);

pop2 = auxdata.PopInterp(lonpop,latpop);

popCost2 = pop2.*AltCost2; % for flights which go over large amounts of

figure()
hold on
plot(t_HeatMin,popCost1)
plot([t_TOSCA_ascent' t_SLEG'+t_TOSCA_ascent(end)],popCost2)
