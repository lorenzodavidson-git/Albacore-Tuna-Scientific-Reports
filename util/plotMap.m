function plotMap(lon,lat,matrix,title1,colbarlabel,colmap,varargin)
    % Check for optional zlimits input
    if nargin > 6 % change to 7 if add landmask
        Zlimits = varargin{1};
    else
        Zlimits = []; % No zlimits provided
    end

    % % backgroundMatrix = double(landmask);  
    % % backgroundMatrix(landmask == 0) = NaN;
    % landmask.land(isnan(landmask.land)) = 0;
    % pcolor(landmask.lon, landmask.lat, landmask.land);
    % cmap = [0.7 0.7 0.7;  % Light blue for ocean
    %     0.7 0.7 0.7]; % Gray for land
    % colormap(gca, cmap); % Light blue for ocean, Gray for land
    % shading flat;
    % clim([0 1]); % Ensures that 1 is mapped to the second color (gray)
    % colorbar('off');
    % hold on
    
    pcolor(lon, lat, matrix);
    shading interp;
    set(gca, 'ylim', [15 62])
    set(gca, 'xlim', [-250 -100])
    if ~isempty(Zlimits)
        clim(Zlimits); % Set color limits if zlimits are given
    end
    grid on
    colorbar;  % Add colorbar
    colBar = colorbar;
    ylabel(colBar, colbarlabel)
    set(get(colBar, 'Label'), 'FontSize', 14);
    xlabel('Longitude');
    ylabel('Latitude');
    title(title1)
    colormap(gca, colmap)
    hold on
    world_coast('linewidth', 2, 'color', 'k')
    shading interp
    set(gca,'Fontsize',14)