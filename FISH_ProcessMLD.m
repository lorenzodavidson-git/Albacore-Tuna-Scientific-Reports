function sst = FISH_ProcessMLD(sst, mld)
% FISH_ProcessMLD - Processes Mixed Layer Depth (MLD) data and interpolates
% it to the Sea Surface Temperature (SST) grid. It then computes the gradient
% angles of the interpolated MLD and adjusts these angles to be within the range [0, 360] degrees.
%
% Syntax: sst = FISH_ProcessMLD(sst, mld)
%
% Inputs:
%   sst - Structure containing the SST data and the corresponding lon/lat grids.
%         This structure must include the following fields:
%           sst.lon - Longitude grid for SST
%           sst.lat - Latitude grid for SST
%
%   mld - Structure containing the MLD data and the corresponding lon/lat grids.
%         This structure must include the following fields:
%           mld.lon - Longitude grid for MLD
%           mld.lat - Latitude grid for MLD
%           mld.seas - 3D array containing MLD data for each month
%
% Outputs:
%   sst - Updated structure with additional fields:
%         sst.mld - Interpolated MLD data on the SST grid
%         sst.mld_angle - Gradient angles of the interpolated MLD data
%         sst.mld_angle_rev - Adjusted gradient angles within the range [0, 360] degrees

% Interpolate MLD to the SST lon, lat
for imon = 1:12
    sst.mld(:,:,imon) = interp2(mld.lon', mld.lat', mld.seas(:,:,imon)', sst.lon, sst.lat);
end

% Compute gradient as a function of month
for imon = 1:12
    [~, ~, sst.mld_angle(:,:,imon), ~] = computeGradient(sst.mld(:,:,imon));
    sst.mld_angle(:,:,imon) = inpaint_nans(sst.mld_angle(:,:,imon), 4);
end

% Reverse the direction by subtracting 180 degrees
sst.mld_angle_rev = sst.mld_angle - 180;

% Find indices where the reversed angles are less than 0 degrees
in_less_than_0 = find(sst.mld_angle_rev < 0);

% Adjust these angles to be within the range [0, 360] by adding 360 degrees
sst.mld_angle_rev(in_less_than_0) = sst.mld_angle_rev(in_less_than_0) + 360;

% Find indices where the reversed angles are greater than or equal to 360 degrees
in_greater_than_360 = find(sst.mld_angle_rev >= 360);

% Adjust these angles to be within the range [0, 360] by subtracting 360 degrees
sst.mld_angle_rev(in_greater_than_360) = sst.mld_angle_rev(in_greater_than_360) - 360;
end
