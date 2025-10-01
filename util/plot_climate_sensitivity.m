function plot_climate_sensitivity(iseas,imon,struct,fish,maxt,mint,name)
    
    ind = imon(:, iseas); ind = ind(ind ~= 0);
    sst = struct.seas;
    mld = struct.mld;
    lon = struct.lon;
    lat = struct.lat;
    
    % Process SST data for the current season
    sstmean = mean(sst(:, :, ind), 3);
    % sstmean2 = interp2(sst.lon', sst.lat', sstmean', mld.lon, mld.lat);
    % sstmean2(sstmean2 > maxt) = nan;
    % sstmean2(sstmean2 < mint) = nan;
    % sstmean2(~isnan(sstmean2)) = 1;
    
    % Process MLD data for the current seaso
    mldmean = mean(mld(:, :, ind), 3);
    mldplot = mldmean;
    mldplot(sstmean > maxt) = nan;
    mldplot(sstmean < mint) = nan;

    p = [-0.020833, 3.3333, 8.5417, 4.3333];
    set(gcf, 'paperposition', p);
    
    % Plot 2D map with SST and MLD data
    FISH_2dmap(mldplot,struct)
    red_blue_colormap
    set(gca, 'ylim', [15 62])
    set(gca, 'xlim', [-250 -100])
    caxis([0 120])
    colorbar;  % Add colorbar
    colBar = colorbar;
    ylabel(colBar, 'MLD (m)')
    colormap(parula);

    % Plot tracks
    in = [];
    for i = 1:6
        ind = find(fish.MON == imon(i, iseas));
        in = [in; ind];
    end
    plot(fish.X(in), fish.Y(in), '.m')
    plot(mean(fish.X(in)),mean(fish.Y(in)),'kx', 'MarkerSize', 14, 'LineWidth', 3)

    % Plot SST contours
    contour(lon, lat, sstmean, [mint maxt], 'linewidth', 6, 'color', 'b');
    title(name)