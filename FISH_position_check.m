function [temp, myflag, theta, p, food] = FISH_position_check(lon1, lat1, imon, sst, flag)
% FISH_position_check: Checks the position of a fish and returns temperature, flag, and direction
%
% This function finds the nearest grid cell in the SST matrix based on the 
% given longitude and latitude, checks the corresponding temperature and probability, 
% and determines if the fish should continue or turn around based on the probability.
%
% Inputs:
%   lon1 - Longitude of the fish's position
%   lat1 - Latitude of the fish's position
%   imon - Month index for the SST data
%   sst  - Structure containing sea surface temperature data and related fields
%   flag - Indicator to choose between land probability and sea probability
%
% Outputs:
%   temp   - Temperature at the fish's position
%   myflag - Boolean indicating if the fish is in favorable conditions (true) or not (false)
%   theta  - Angle of the temperature gradient direction (in radians)
%   p      - Probability value at the fish's position

% Find the nearest grid cell in the SST matrix
[~, lonIndex] = min(abs(sst.lon(:, 1) - lon1));
[~, latIndex] = min(abs(sst.lat(1, :) - lat1));

% Get the temperature and gradient angle for the specified month
temp = sst.seas(lonIndex, latIndex, imon);
gangle = sst.grad_angle(lonIndex, latIndex, imon);
food.prod = sst.lme_prod(lonIndex, latIndex, imon);
food.theta = deg2rad(sst.lme_angle(lonIndex, latIndex));


food.mld = sst.mld(lonIndex, latIndex, imon);
food.mld_theta = deg2rad(sst.mld_angle(lonIndex, latIndex,1));
food.mld_theta_rev = deg2rad(sst.mld_angle_rev(lonIndex, latIndex,imon));

% Determine the probability based on the flag
if flag == 0
    p = sst.prob_land(lonIndex, latIndex);
else
    p = sst.prob(lonIndex, latIndex, imon);
end

% Determine if the fish is in favorable conditions (myflag)
% myflag is true if the fish is in favorable conditions and false otherwise
myflag = rand < abs(p) / 100;

% Convert the gradient angle to radians
theta = deg2rad(gangle);

% Additional commented-out logic for further checks on temperature
% Uncomment and adjust as needed based on specific requirements
% myflag = 1;
% if temp < 10.8089
%     myflag = false;
% end
% if temp > 23.3153
%     myflag = false;
% end
% if isnan(temp)
%     gangle = 180;
%     myflag = false;
% end

% Example logic for calculating components of a movement vector
% Uncomment and adjust as needed based on specific requirements
% r = 1; % Example magnitude
% dx = r * cos(theta) * sign(p);
% dy = r * sin(theta) * sign(p);

end
