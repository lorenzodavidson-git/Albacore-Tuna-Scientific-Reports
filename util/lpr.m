function lpr(varargin)
% lpr: Set properties for all axes in the current figure and print to an EPS file
%
% This function adjusts the properties of all axes in the current figure, 
% setting the line width and layer properties, and then prints the figure 
% to an EPS file with high resolution.
%
% Inputs:
%   filename - Optional argument specifying the filename for the EPS file.
%              If not provided, the default filename is 'tmp.eps'.

% Check if a filename is provided as an input argument
if nargin > 0
    filename = varargin{1};
else
    filename = 'tmp.eps'; % Default filename
end

% Find all axes in the current figure
allAxesInFigure = findall(gcf, 'type', 'axes');

% Set properties for each axis
for i = 1:length(allAxesInFigure)
    axish = allAxesInFigure(i);
    set(axish, 'LineWidth', 1); % Set the line width
    set(axish, 'Layer', 'top'); % Set the layer to 'top'
end

% Print the figure to an EPS file with specified options
print('-depsc2', '-painters', '-loose', '-r1200', filename);

end
