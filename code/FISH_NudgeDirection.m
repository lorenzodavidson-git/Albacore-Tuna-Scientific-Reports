function [nudgedU, nudgedV] = FISH_NudgeDirection(u, v, desiredAngle, weight)
% FISH_NudgeDirection - Nudge the direction of velocity vectors towards a desired angle
%
% Syntax: [nudgedU, nudgedV] = FISH_NudgeDirection(u, v, desiredAngle, weight)
%
% Inputs:
%   u           - x-component of the original velocity vector
%   v           - y-component of the original velocity vector
%   desiredAngle - Desired direction to nudge the vectors towards (in degrees)
%   weight      - Weighting factor (0 to 1, where 0 means no nudging and 1 means fully nudged)
%
% Outputs:
%   nudgedU     - x-component of the nudged velocity vector
%   nudgedV     - y-component of the nudged velocity vector

% Calculate the magnitude of the original velocity vector
amp = sqrt(u^2 + v^2);

% Convert the desired angle to radians
desiredAngleRad = deg2rad(desiredAngle);

% Calculate the desired direction vector components
desiredU = cos(desiredAngleRad);
desiredV = sin(desiredAngleRad);

% Nudge the random vectors towards the desired direction
nudgedU = (1 - weight) * u + weight * desiredU;
nudgedV = (1 - weight) * v + weight * desiredV;

% Normalize the nudged vectors to maintain the original magnitude
normFactor = sqrt(nudgedU^2 + nudgedV^2);
nudgedU = nudgedU / normFactor * amp;
nudgedV = nudgedV / normFactor * amp;

end


