function [sst] = FISH_ComputeProbability_Cont(sst, fish)
% FISH_ComputeProbability: Computes probability maps based on surface temperature data
%
% This function processes surface temperature data, computes temperature gradients, 
% and calculates probability maps based on the temperature data and fish tracking data.
%
% Inputs:
%   sst  - Structure containing sea surface temperature data
%   fish - Structure containing fish tracking data
%
% Outputs:
%   sst  - Updated structure with computed probability maps and gradients

% Flag for plotting (0 = no plotting, 1 = plotting)
PLOTTING = 0;

% Process surface temperature data for each month
% for imon = 1:12
%     in = find(imon == sst.month);
%     sst.seas(:,:,imon) = mean(sst.data(:,:,in), 3);
% end
sst.mean = mean(sst.seas, 3);

% Compute temperature gradient, direction, and amplitude
[Gx, Gy, grad_angle, amplitude] = computeGradient(sst.mean);

% Define indices for quiver plot subsampling
[I, J] = size(sst.lon);
dx = 40; I = 1:dx:I; J = 1:dx:J;

% Plot mean surface temperature and temperature gradient direction if plotting is enabled
if PLOTTING == 1
    figure(1)
    clf
    subplot(2,1,1)
    plot_2dmap(sst.mean, sst)
    title('Mean Surface Temperature from 1980-2024')
    set(gca, 'FontSize', 14)
    
    subplot(2,1,2)
    plot_2dmap(grad_angle, sst)
    quiver(sst.lon(I,J), sst.lat(I,J), Gx(I,J), Gy(I,J), 4, 'k')
    colormap("parula")
    title('Direction of Temperature Gradient (angle 90 degree = North)')
    set(gca, 'FontSize', 14)
end

% Process gradient for each month
for imon = 1:12
    [Gx(:,:,imon), Gy(:,:,imon), sst.grad_angle(:,:,imon), ~] = ...
        computeGradient(sst.seas(:,:,imon));
end


% Initialize temperature and probability matrices
TEMP = sst.seas;
PROB = sst.seas;
PROB(:) = 0;

% Compute mean and standard deviation of fish tracking data
mu = meanNaN(fish.ta.data, 1)+0.0;
sigma = stdNaN(fish.ta.data, 1);

disp ('Using new SST probability')

%% Assign probability values based on temperature thresholds

% within 1 STD give a probability of 100%
in = find(TEMP(:) < mu + sigma * 1 & TEMP(:) > mu - sigma * 1); %%%%%%% replaced stdNaN(fish.ta.data, 1) with sigma
PROB(in) = 100;

in = find(TEMP(:) >= mu + sigma * 1); PROB(in) = -85;
in = find(TEMP(:) >= mu + sigma * 2); PROB(in) = -65;
in = find(TEMP(:) >= mu + sigma * 3); PROB(in) = -20;
in = find(TEMP(:) >= mu + sigma * 4); PROB(in) = -5;

in = find(TEMP <= mu - sigma * 1); PROB(in) = 85; %%%%%%% add (:) ?
in = find(TEMP <= mu - sigma * 2); PROB(in) = 65; %%%%%%% add (:) ?
in = find(TEMP <= mu - sigma * 3); PROB(in) = 20; %%%%%%% add (:) ?
in = find(TEMP(:) <= mu - sigma * 4); PROB(in) = 5;

% Update sst structure with computed probabilities
sst.prob = PROB;

% Handle NaN values in seas and grad_angle
in = find(isnan(sst.seas));
sst.prob(in) = 1;

[I,J,N] = size(sst.grad_angle);
for n = 1 : N
    tmp = sst.grad_angle(:,:,n);
    in = find(isnan(tmp));
    tmp(in) = sst.mask_grad(in);
    sst.grad_angle(:,:,n) = tmp;
end


% Compute land mask probability
prob = sst.mask;
prob(~isnan(prob)) = 100;
prob(isnan(prob)) = 1;
sst.prob_land = prob;

end

