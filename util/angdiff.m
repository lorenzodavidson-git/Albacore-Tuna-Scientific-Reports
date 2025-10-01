function diff = angdiff(theta2, theta1)
    % Compute minimal difference between two angles in radians
    diff = mod(theta2 - theta1 + pi, 2*pi) - pi;
end