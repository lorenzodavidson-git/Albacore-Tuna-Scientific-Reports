function [Gx, Gy, radians, amplitude] = computeGradient(A)

[rows, cols] = size(A);

% Initialize gradient matrices
Gx = zeros(rows, cols);
Gy = zeros(rows, cols);

% Compute gradient in x-direction (rows)
Gx(2:rows-1, :) = (A(3:rows, :) - A(1:rows-2, :)) / 2;
% For the first and last row, use forward and backward difference
Gx(1, :) = A(2, :) - A(1, :);
Gx(rows, :) = A(rows, :) - A(rows-1, :);

% Compute gradient in y-direction (columns)
Gy(:, 2:cols-1) = (A(:, 3:cols) - A(:, 1:cols-2)) / 2;
% For the first and last column, use forward and backward difference
Gy(:, 1) = A(:, 2) - A(:, 1);
Gy(:, cols) = A(:, cols) - A(:, cols-1);
radians = atan2d(Gy, Gx); % atan2d is used to handle all quadrants
in = find(radians < 0);
radians(in) = radians(in) + 360;
amplitude = sqrt(Gx.^2 + Gy.^2);


