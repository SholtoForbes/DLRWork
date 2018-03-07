% World population interpolator
% Sholto Forbes-Spyratos 2018
clear all

% [popdata,tiffdata] = geotiffread('pop_lowres');
% [popdata,tiffdata] = geotiffread('gpw_v4_population_count_rev10_2020_1_deg');
[popdata,tiffdata] = geotiffread('gpw_v4_population_count_rev10_2020_2pt5_min');
popdata(popdata < 0) = 0;

popdata = flipud(popdata);

LonList = tiffdata.LongitudeLimits(1):(tiffdata.LongitudeLimits(end)-tiffdata.LongitudeLimits(1))/(tiffdata.RasterSize(2)-1):tiffdata.LongitudeLimits(end);

LatList = tiffdata.LatitudeLimits(1):(tiffdata.LatitudeLimits(end)-tiffdata.LatitudeLimits(1))/(tiffdata.RasterSize(1)-1):tiffdata.LatitudeLimits(end);

[LonGrid,LatGrid] = meshgrid(LonList,LatList);

% PopInterp = griddedInterpolant(LonGrid',LatGrid',popdata','spline');
PopInterp = griddedInterpolant(LonGrid',LatGrid',popdata','linear');


figure()
contour(LonGrid,LatGrid,popdata,100)

figure()
contour(LonGrid,LatGrid,PopInterp(LonGrid,LatGrid),100)

save('PopInterp','PopInterp');