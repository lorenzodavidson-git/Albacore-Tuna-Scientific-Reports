function make_box(xlim, ylim, varargin)
% Make_box: Draw a rectangular box on the current plot
%
% This function draws a rectangular box on the current plot using the specified
% x and y limits. The box can be customized with an optional color argument.
%
% Inputs:
%   xlim     - A vector specifying the x-axis limits [xmin xmax].
%   ylim     - A vector specifying the y-axis limits [ymin ymax].
%   varargin - Optional argument specifying the color of the box.
%              If not provided, the default color is blue ('b').

% Define the corners of the box based on the provided limits
l1 = [min(xlim) min(xlim) max(xlim) max(xlim) min(xlim)];
l2 = [min(ylim) max(ylim) max(ylim) min(ylim) min(ylim)];

% Hold the current plot to add the box
hold on

% Check if a color is provided as an input argument
if nargin > 2
    c = varargin{1};    
    plot(l1, l2, 'linewidth', 4, 'color', c); % Plot the box with the specified color
else
    plot(l1, l2, 'linewidth', 4, 'color', 'k'); % Plot the box with the default color (black)
end

end
