function [sim_struct] = FISH_run_simulation(projection,fish2d,fish,type)

% Initialize figure
fig7 = figure;
set(gcf, 'paperposition', [-0.020833, 3.3333, 6.5417, 4.3333]);
clf

% Plot land map
FISH_2dmap(projection.land, projection)
set(gca, 'ylim', [15 62])
colormap([0.8, 0.8, 0.8])
colorbar off

% Initialize data storage
MON = [];
U = [];
V = [];
X = [];
Y = [];
TEMP = [];
MLD = [];

for n = 1:10
    for i = 1:12
        % Extract fish data
        lon = fish2d(i).lon(1);
        lat = fish2d(i).lat(1);
        u = fish2d(i).u;
        v = fish2d(i).v;
        dt = fish2d(i).time(2:end) - fish2d(i).time(1:end-1);
        time = (fish2d(i).time(2:end) + fish2d(i).time(1:end-1)) / 2;
        month = str2num(datestr(time, 'mm'));

        % Simulate fish swimming
        [x, y, tf, u, v, month, mld_fish] = FISH_simulate_all_edit3(lon, lat, u', v', dt, time, projection, fish, type);
        
        MON = [MON; month];

        % Store results
        U = [U; u];
        V = [V; v];
        TEMP = [TEMP; tf(1:end-1)];
        MLD = [MLD; mld_fish(1:end)];
        X = [X; x(1:end-1)];
        Y = [Y; y(1:end-1)];
        
        % Plot fish trajectory
        plot(x, y, 'LineWidth', 1)
    end
end

% Plot fish trajectories in black
for i = 1:12
    plot(fish2d(i).lon, fish2d(i).lat, 'LineWidth', 2, 'color', 'k')
end

sim_struct = struct();
sim_struct.MON = MON;
sim_struct.U = U;
sim_struct.V = V;
sim_struct.TEMP = TEMP;
sim_struct.MLD = MLD;
sim_struct.X = X;
sim_struct.Y = Y;