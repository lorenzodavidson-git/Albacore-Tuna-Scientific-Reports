function FISH_2dmap (field2d, data)
% plot_2dmap (field2d, data)
% June 27, 2023 by Connie Erdozain
% field2d: is a 2D array contianing the spatial field to plot
% data: is structured variable containing the lon, lat information
% data = 
%      lon: [180×89 double]
%      lat: [180×89 double]


%clf;
pcolor(data.lon, data.lat, field2d);
colorbar;
hold on
world_coast('linewidth', 2, 'color', 'k')
shading interp
set(gca,'Fontsize',14)