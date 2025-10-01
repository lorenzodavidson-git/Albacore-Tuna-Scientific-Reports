function lpr_png(varargin)
% lpr_png: Set properties for all axes in the current figure and print to a PNG file
%
% This function adjusts the properties of all axes in the current figure,
% setting the line width and layer properties, and then prints the figure
% to a PNG file with high resolution.
%
% Inputs:
%   filename - Optional argument specifying the filename for the PNG file.
%              If not provided, the default filename is 'tmp.png'.

% Check if a filename is provided as an input argument
if nargin > 0
    filename = varargin{1};
else
    filename = 'tmp.png'; % Default filename
end

% Find all axes in the current figure
allAxesInFigure = findall(gcf, 'type', 'axes');

% Set properties for each axis
for i = 1:length(allAxesInFigure)
    axish = allAxesInFigure(i);
    set(axish, 'LineWidth', 1); % Set the line width
    set(axish, 'Layer', 'top'); % Set the layer to 'top'
end

% Print the figure to a PNG file with 300 dpi resolution
% name = fullfile('/Desktop','Post-Brown Research','Fish Paper','figs',filename);
print(filename, '-dpng', '-r300');

end
