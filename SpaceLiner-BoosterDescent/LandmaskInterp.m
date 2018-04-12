% function [LandSpline] = LandmaskInterp()

lon = [-180:1:180];
lat = [-90:1:90];

[lonGrid,latGrid] = meshgrid(lon,lat);

for i = 1:length(lon)
    i
    for j = 1:length(lat)
        LandGrid(j,i) = landmask(latGrid(j,i),lonGrid(j,i),90);
    end
end

LandSpline = griddedInterpolant(lonGrid',latGrid',double(LandGrid)','nearest');

save('LandSpline','LandSpline')

% end

