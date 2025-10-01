function sst = FISH_CCS_LME_Region(sst, prod)
% FISH_CCS_LME_Region This function modifies the input SST data by applying
% a mask for the California Current System (CCS) Large Marine Ecosystem (LME) region.
%
% Inputs:
%   sst - a struct containing the sea surface temperature (SST) data and
%         associated mask, latitude, and longitude grids.
%         Fields required in the sst struct:
%         - sst.mask: a mask indicating valid SST data points
%         - sst.lat: latitude grid
%         - sst.lon: longitude grid
%   prod - a 1x12 array containing monthly productivity data to be added to
%          the SST data.
%
% Outputs:
%   sst - the modified input struct with additional fields:
%         - sst.lme_mask: mask for the CCS LME region
%         - sst.lme_prod: modified SST data with applied mask and added productivity data
%         - sst.lme_angle: angle to the nearest blue areas

% Define region of CCS LME using longitude and latitude coordinates
x = [-116.275870483041, -135.878752109265, -136.850795826268, ...
     -137.984846829438, -139.280905118775, -139.60491969111, ...
     -139.928934263444, -138.794883260274, -137.012803112435, ...
     -135.554737536931, -132.638606385922, -130.208497093415, ...
     -127.13035865624, -124.700249363732, -120.650067209554, ...
     -115.951855910707, -112.387695615029, -106.879447885346, ...
     -108.499520747018, -111.253644611859, -116.275870483041];

y = [55.9099164341254, 55.9099164341254, 53.3059914690367, ...
     51.1026703447309, 48.9995001806209, 47.1967828970979, ...
     44.4927069718135, 41.0875743251591, 38.4836493600705, ...
     36.7810830367433, 35.078516713416, 32.975346549306, ...
     30.8721763851959, 29.6703648628473, 28.5687043006944, ...
     27.6673456589329, 27.3668927783458, 28.0679494997158, ...
     41.3880272057463, 49.5002549815994, 55.9099164341254];

% Initialize a new mask with NaN values
newMask = sst.mask * nan;

% Extract latitude and longitude grids from the sst struct
latGrid = sst.lat;
lonGrid = sst.lon;

% Iterate over each point in the mask
for i = 1:size(newMask, 1)
    for j = 1:size(newMask, 2)
        % Check if the point is inside the polygon defined by x and y
        if inpolygon(lonGrid(i, j), latGrid(i, j), x, y)
            % Set the mask value to 1 for points inside the polygon
            newMask(i, j) = 1;
        end
    end
end

% Update the sst struct with the new LME mask

sst.lme_mask = newMask .* sst.mask;

% Apply the productivity data to the SST data for each month
for imon = 1:12
    sst.lme_prod(:, :, imon) = sst.lme_mask + prod(imon) - 1;
    
end

sst.lme_prod(isnan(sst.lme_prod)) = 0;
sst.lme_prod = sst.lme_prod.*repmat(sst.mask,[1 1 12]);

% Optionally, you can visualize the data using a custom mapping function
% mfig
% FISH_2dmap(sst.lme_prod(:, :, 1), sst)

%% Find angle of swim direction towards CCS

% Convert the image to double for processing
mapData = sst.lme_mask;
mapData(isnan(mapData)) = 0;

% Identify the blue areas (adjust the threshold according to your map)
% Assuming blue areas are marked with specific values
blueThreshold = 1; % Adjust this value according to your data
blueMask = mapData == blueThreshold;

% Get the size of the map
[rows, cols] = size(mapData);

% Compute the distance transform and nearest neighbor indices
[dist, nearestIdx] = bwdist(blueMask);

% Initialize array to store the nearest direction for each point
nearestDirection = nan(rows, cols);

% Get the coordinates of all points
[X, Y] = meshgrid(1:cols, 1:rows);

% Get the coordinates of the nearest blue areas
[nearestBlueRows, nearestBlueCols] = ind2sub(size(mapData), nearestIdx);

% Calculate the direction to the nearest blue area
deltaY = nearestBlueRows - Y;
deltaX = nearestBlueCols - X;
angles = atan2d(deltaY, deltaX);

% Adjust angles to be in the range [0, 360]
angles(angles < 0) = angles(angles < 0) + 360;

% Set directions only for non-blue areas
nearestDirection(~blueMask) = angles(~blueMask);

sst.lme_angle = nearestDirection;





end
