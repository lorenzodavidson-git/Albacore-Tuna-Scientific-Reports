function [circ_median_present_days, circ_median_future_days, circ_diff_days, ci] = circular_bootstrap_diff(present_data, future_data, nBoot)
% circular_bootstrap_diff
% Computes circular medians and bootstrapped confidence interval of the difference.
%
% Inputs:
%   present_data - vector of present scenario day-of-year data
%   future_data  - vector of future scenario day-of-year data
%   nBoot        - number of bootstrap iterations (default = 5000)
%
% Outputs:
%   circ_median_present_days - circular median of present data (days)
%   circ_median_future_days - circular median of future data (days)
%   circ_diff_days          - difference in circular medians (days)
%   ci                      - 95% confidence interval for difference (days)

if nargin < 3
    nBoot = 5000;
end

% Convert to radians
present_rad = present_data / 365 * 2 * pi;
future_rad = future_data / 365 * 2 * pi;

% Compute circular mean (proxy for median)
circ_mean_present = atan2(mean(sin(present_rad)), mean(cos(present_rad)));
circ_mean_future = atan2(mean(sin(future_rad)), mean(cos(future_rad)));

% Convert to days
circ_median_present_days = mod(circ_mean_present, 2*pi) / (2*pi) * 365;
circ_median_future_days = mod(circ_mean_future, 2*pi) / (2*pi) * 365;

% Inline angdiff function
angdiff = @(theta2, theta1) mod(theta2 - theta1 + pi, 2*pi) - pi;

% Compute circular difference in days
circ_diff_days = angdiff(circ_mean_future, circ_mean_present) / (2*pi) * 365;

% Bootstrap
boot_diff_days = zeros(nBoot,1);
n_present = numel(present_data);
n_future = numel(future_data);

for i = 1:nBoot
    % Resample with replacement
    resample_present = present_data(randi(n_present, n_present, 1));
    resample_future = future_data(randi(n_future, n_future, 1));
    
    % Convert to radians
    resample_present_rad = resample_present / 365 * 2 * pi;
    resample_future_rad = resample_future / 365 * 2 * pi;
    
    % Compute circular mean (proxy for median)
    mean_present = atan2(mean(sin(resample_present_rad)), mean(cos(resample_present_rad)));
    mean_future = atan2(mean(sin(resample_future_rad)), mean(cos(resample_future_rad)));
    
    % Difference (circular), convert to days
    diff_rad = angdiff(mean_future, mean_present);
    boot_diff_days(i) = diff_rad / (2*pi) * 365;
end

% Compute 95% confidence interval
ci = prctile(boot_diff_days, [2.5 97.5]);

% Print output
fprintf('Circular median (present): %.2f days\n', circ_median_present_days);
fprintf('Circular median (future): %.2f days\n', circ_median_future_days);
fprintf('Circular difference: %.2f days\n', circ_diff_days);
fprintf('Bootstrap CI for difference: %.2f to %.2f days\n', ci(1), ci(2));

end