function [shifted_doy, mean_doy] = circular_doy(doy_vector)

% Convert to radians and find mean doy
radians = doy_vector / 365 * 2 * pi;
x = cos(radians);
y = sin(radians);
mean_angle = atan2(mean(y,'omitnan'), mean(x,'omitnan'));
mean_doy = mod(mean_angle, 2*pi) / (2*pi) * 365;

shifted_doy = mod(doy_vector - mean_doy + 182.5, 365);
end
