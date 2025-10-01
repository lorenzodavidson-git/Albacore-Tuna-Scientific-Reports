%% Add all data and paths

addpath util
addpath data
addpath fish_data
addpath code


% Climate

sst = Load_GLORYS_sst;
mask = sst.mean;
in=find(isnan(mask));
mask(~isnan(mask))=1;
sst.mask = mask;
mask(:)=nan; mask(in)=1;
sst.land = mask;

load GLORYS_mld_seasonal_smooth.mat

prod =[-1.2 -1.2 -1.1 -.5 0 0.6 1.2 1.4 1.3 0.7 -.1 -.85];
prod = prod-mean(prod); prod= prod/std(prod);

% Fish stats for SST gradient processing

load albacore_fish.mat % fish --> fish_all
fish_all = fish;
clear fish

load albacore_fish2d % fish2d

load fish_stats.mat % fish

%% Process Environmental Data

% SST
sst = FISH_ComputeMaskGradient(sst);
sst = FISH_CCS_LME_Region (sst, prod);
sst = FISH_ProcessMLD(sst,mld);
sst = FISH_ComputeProbability_Cont(sst, fish);

% MLD
lat_max = max(sst.lat(1,:)); lat_min = min(sst.lat(1,:));
lon_max = max(sst.lon(:,1)); lon_min = min(sst.lon(:,1));

[mld.lon,mld.lat,mld.seas] = clipData(mld.seas,mld.lon(:,1),mld.lat(1,:),lon_min,lon_max,lat_min,lat_max);
[~,~,mld.mask] = clipData(mld.mask,mld.lon(:,1),mld.lat(1,:),lon_min,lon_max,lat_min,lat_max);
[~,~,mld.mean] = clipData(mld.mean,mld.lon(:,1),mld.lat(1,:),lon_min,lon_max,lat_min,lat_max);
mld.land = mask;
for i = 1:12
    layer = mld.seas(:,:,i);
    layer(mld.land == 1) = NaN;
    mld.seas(:,:,i) = layer;
end
mld.mean(sst.land == 1) = NaN;
for i = 1:12
    layer = sst.mld(:,:,i);
    layer(sst.land == 1) = NaN;
    sst.mld(:,:,i) = layer;
end

%% Process Fish Data

% Relevant data for migratory fish for 2D spatial analysis
igood = [1 2 3 5 9 11 13 14 15 17 21 22];
ibad = [4 6 7 8 10 12 16 18 19 20];
fish2d = fish2d(igood);

% All relevant data for migratory fish for vertical analysis
tags_12 = [394, 396, 1045, 1464, 1987, 2082, 2381, 2393, 2398, 2942, 1090251, 1090269];
in = find(ismember(fish_all.tags, tags_12));
mig_fish = struct();
mig_fish.tag = fish_all.tags(in);
mig_fish.lon = fish_all.lon(in);
mig_fish.lat = fish_all.lat(in);
mig_fish.year = fish_all.year(in);
mig_fish.month = fish_all.month(in);
mig_fish.day = fish_all.day(in);
mig_fish.time = fish_all.datenum(in);
mig_fish.depth = fish_all.depth(in);
mig_fish.dn = fish_all.dn(in);

for t=1:length(mig_fish.depth)
    imon = mig_fish.month(t);
    [~, lonIndex] = min(abs(mld.lon(:, 1) - mig_fish.lon(t)));
    [~, latIndex] = min(abs(mld.lat(1, :) - mig_fish.lat(t)));
    mig_fish.mld(t) = mld.seas(lonIndex,latIndex,imon);
end
mig_fish.mld = mig_fish.mld(:);

clear i ibad igood in layer imon latIndex lonIndex t

%% Calculate daily, daytime and nighttime depths & MLD for mig_fish

num_days = 365;
bin_daymonth = datetime(2023, 1, 1):datetime(2023, 12, 31);

for i = 1:num_days
    % Extract month and day
    daymonth = bin_daymonth(i);
    imonth = month(daymonth);
    iday = day(daymonth);

    % Find indices for the current day-month
    in = (mig_fish.month == imonth) & (mig_fish.day == iday);

    % Separate into day and night
    in_day = strcmp(mig_fish.dn(in), 'day');
    in_night = strcmp(mig_fish.dn(in), 'night');

    % Get corresponding data
    depth = mig_fish.depth(in);
    mld_data = mig_fish.mld(in);
    depth_day = depth(in_day);
    depth_night = depth(in_night);

    % Calculate metrics
    mig_fish.daily.lon(i) = mean(mig_fish.lon(in), 'omitnan');
    mig_fish.daily.lat(i) = mean(mig_fish.lat(in), 'omitnan');

    mig_fish.daily.depth(i) = mean(depth, 'omitnan');
    mig_fish.daily.mld(i) = mean(mld_data, 'omitnan');
    mig_fish.daily.day_depth(i) = mean(depth_day, 'omitnan');
    mig_fish.daily.night_depth(i) = mean(depth_night, 'omitnan');

    mig_fish.daily.depth_median(i) = median(depth, 'omitnan');
    mig_fish.daily.day_depth_median(i) = median(depth_day, 'omitnan');
    mig_fish.daily.night_depth_median(i) = median(depth_night, 'omitnan');

    mig_fish.daily.mld_median(i) = median(mld_data, 'omitnan');
    mig_fish.daily.day_mld_median(i) = median(mld_data(in_day), 'omitnan');
    mig_fish.daily.night_mld_median(i) = median(mld_data(in_night), 'omitnan');

    mig_fish.daily.depth_90(i) = prctile(depth, 90);
    mig_fish.daily.day_depth_90(i) = prctile(depth_day, 90);
    mig_fish.daily.night_depth_90(i) = prctile(depth_night, 90);

    mig_fish.daily.depth_top10(i) = mean(depth(depth >= prctile(depth, 90)));
    mig_fish.daily.day_depth_top10(i) = mean(depth_day(depth_day >= prctile(depth_day, 90)));
    mig_fish.daily.night_depth_top10(i) = mean(depth_night(depth_night >= prctile(depth_night, 90)));
end

for imon = 1:12
    in = mig_fish.month == imon; % Find indices for the current month
    in_day = strcmp(mig_fish.dn(in), 'day');
    in_night = strcmp(mig_fish.dn(in), 'night');

    depth = mig_fish.depth(in);
    mld_data = mig_fish.mld(in);

    mig_fish.seas.lon(imon) = mean(depth, 'omitnan'); % Average depth
    mig_fish.seas.lat(imon) = mean(depth, 'omitnan'); % Average depth
    mig_fish.seas.depth(imon) = mean(depth, 'omitnan'); % Average depth
    mig_fish.seas.mld(imon) = mean(mld_data, 'omitnan'); % Average depth
    mig_fish.seas.day_depth(imon) = mean(depth(in_day), 'omitnan'); % Average depth
    mig_fish.seas.night_depth(imon) = mean(depth(in_night), 'omitnan'); % Average depth
    mig_fish.seas.depth_median(imon) = median(depth, 'omitnan');
    mig_fish.seas.day_depth_median(imon) = median(depth(in_day), 'omitnan');
    mig_fish.seas.night_depth_median(imon) = median(depth(in_night), 'omitnan');
end

%% Create Indiv_fish for analyzing seasonality of individual tags (columns 1-12), average (column 13) and std (column 14)

num_days = 365;
bin_daymonth = datetime(2023, 1, 1):datetime(2023, 12, 31);
tags_12 = [394, 396, 1045, 1464, 1987, 2082, 2381, 2393, 2398, 2942, 1090251, 1090269];

for i = 1:num_days
    % Extract month and day
    daymonth = bin_daymonth(i);
    imonth = month(daymonth);
    iday = day(daymonth);

    for j = 1:length(tags_12)
        in = (mig_fish.tag == tags_12(j)) & (mig_fish.month == imonth) & (mig_fish.day == iday);

        lon_in = mig_fish.lon(in);
        lat_in = mig_fish.lat(in);
        mld_in = mig_fish.mld(in);
        depth_in = mig_fish.depth(in);

        in_day = strcmp(mig_fish.dn(in), 'day');
        in_night = strcmp(mig_fish.dn(in), 'night');

        depth_day_in = depth_in(in_day);
        depth_night_in = depth_in(in_night);

        % Create fields
        indiv_fish.lon(i,j) = mean(lon_in, 'omitnan');
        indiv_fish.lat(i,j) = mean(lat_in, 'omitnan');
        indiv_fish.mld(i,j) = mean(mld_in, 'omitnan');
        indiv_fish.depth(i,j) = mean(depth_in, 'omitnan');
        indiv_fish.day_depth(i,j) = mean(depth_day_in, 'omitnan');
        indiv_fish.night_depth(i,j) = mean(depth_night_in, 'omitnan');

        indiv_fish.depth_median(i,j) = prctile(depth_in, 50);
        indiv_fish.day_depth_median(i,j) = prctile(depth_day_in, 50);
        indiv_fish.night_depth_median(i,j) = prctile(depth_night_in, 50);

        indiv_fish.depth_top10(i,j) = mean(depth_in(depth_in >= prctile(depth_in, 90)));
        indiv_fish.day_depth_top10(i,j) = mean(depth_day_in(depth_day_in >= prctile(depth_day_in, 90)));
        indiv_fish.night_depth_top10(i,j) = mean(depth_night_in(depth_night_in >= prctile(depth_night_in, 90)));
    end

    names = fieldnames(indiv_fish);
    for x = 1:length(names)
        current_field = names{x};  % Extract field name properly
        indiv_fish.(current_field)(i,13) = mean(indiv_fish.(current_field)(i,1:length(tags_12)), 'omitnan');
        indiv_fish.(current_field)(i,14) = std(indiv_fish.(current_field)(i,1:length(tags_12)), 'omitnan');
    end

end

%% Create fish_vertical to contain daily averaged depth and MLD for spatial maps

tags_12 = [394, 396, 1045, 1464, 1987, 2082, 2381, 2393, 2398, 2942, 1090251, 1090269];

% Calculate average depth for every day for each tuna
fish_vertical = struct();
for i = 1:length(tags_12)
    subset = struct();
    in = find(ismember(mig_fish.tag, tags_12(i)));
    subset.lon = mig_fish.lon(in);
    subset.lat = mig_fish.lat(in);
    subset.depth = mig_fish.depth(in);
    subset.mld = mig_fish.mld(in);
    subset.time = mig_fish.time(in);
    subset.dn = mig_fish.dn(in);
    
    days = floor(subset.time);
    unique_days = unique(days);
    daytime = strcmp(subset.dn, 'day');   % Logical array for "day"
    nighttime = strcmp(subset.dn, 'night'); % Logical array for "night"

        for j = 1:length(unique_days)
            daily_indices = (days == unique_days(j)); % Indices for the current day

            in_day = strcmp(subset.dn(daily_indices), 'day');
            in_night = strcmp(subset.dn(daily_indices), 'night');

            % Get corresponding data
            lon = subset.lon(daily_indices);
            lat = subset.lat(daily_indices);
            depth_subset = subset.depth(daily_indices);
            mld_subset = subset.mld(daily_indices);
            depth_day = depth_subset(in_day);
            depth_night = depth_subset(in_night);

            subset.daily_lon(j) = mean(lon); subset.daily_lon = subset.daily_lon(:);
            subset.daily_lat(j) = mean(lat); subset.daily_lat = subset.daily_lat(:);
            subset.daily_depth(j) = mean(depth_subset); subset.daily_depth = subset.daily_depth(:);
            subset.daily_mld(j) = mean(mld_subset); subset.daily_mld = subset.daily_mld(:);

            subset.day_depth(j) = mean(depth_day); subset.day_depth = subset.day_depth(:);
            subset.night_depth(j) = mean(depth_night); subset.night_depth = subset.night_depth(:);

            subset.daily_depth_median(j) = prctile(depth_subset, 50); subset.daily_depth_median = subset.daily_depth_median(:);
            subset.day_depth_median(j) = prctile(depth_day, 50); subset.day_depth_median = subset.day_depth_median(:);
            subset.night_depth_median(j) = prctile(depth_night, 50); subset.night_depth_median = subset.night_depth_median(:);

            subset.daily_depth_top10(j) = mean(depth_subset(depth_subset >= prctile(depth_subset, 90))); subset.daily_depth_top10 = subset.daily_depth_top10(:);
            subset.day_depth_top10(j) = mean(depth_day(depth_day >= prctile(depth_day, 90))); subset.day_depth_top10 = subset.day_depth_top10(:);
            subset.night_depth_top10(j) = mean(depth_night(depth_night >= prctile(depth_night, 90))); subset.night_depth_top10 = subset.night_depth_top10(:);
        end
    
    subset.unique_day = unique_days;
    subset.day = datetime(unique_days, 'ConvertFrom', 'datenum');
    
    field_name = sprintf('tag_%d', tags_12(i));
    fish_vertical.(field_name) = subset;
end

fish_vertical.daily = struct('day_of_year', (1:365)', 'daily_lon', nan(365, 13), 'daily_lat', nan(365, 13), 'daily_depth', nan(365, 13), ...
    'daily_mld', nan(365, 13), 'day_depth', nan(365, 13), 'night_depth', nan(365, 13));

%% SST and MLD Projection Anomaly Data

addpath data/Unprocessed_SSTMLD_Projection_Anomalies

% SST projections
fname = 'SST_SSP245_annualAnom_19852014_20702099.nc';
readncfile
lon = lon - 360;
SST_anom = struct();
[SST_anom.lon,SST_anom.lat,SST_anom.year] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);

fname = 'SST_SSP245_annualAnom_19852014_20702099_JFM.nc';
readncfile
lon = lon - 360;
[~,~,SST_anom.jfm] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);
fname = 'SST_SSP245_annualAnom_19852014_20702099_AMJ.nc';
readncfile
lon = lon - 360;
[~,~,SST_anom.amj] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);
fname = 'SST_SSP245_annualAnom_19852014_20702099_JAS.nc';
readncfile
lon = lon - 360;
[~,~,SST_anom.jas] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);
fname = 'SST_SSP245_annualAnom_19852014_20702099_OND.nc';
readncfile
lon = lon - 360;
[~,~,SST_anom.ond] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);

SST_anom.lon_interp = sst.lon; SST_anom.lat_interp = sst.lat; 
SST_anom = interpolateData(SST_anom); % Interpolate

% MLD projections
fname = 'MLD_SSP245_annualAnom_19852014_20702099.nc';
readncfile
lon = lon - 360;
MLD_anom = struct();
[MLD_anom.lon,MLD_anom.lat,MLD_anom.year] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);

fname = 'MLD_SSP245_annualAnom_19852014_20702099_JFM.nc';
readncfile
lon = lon - 360;
[~,~,MLD_anom.jfm] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);
fname = 'MLD_SSP245_annualAnom_19852014_20702099_AMJ.nc';
readncfile
lon = lon - 360;
[~,~,MLD_anom.amj] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);
fname = 'MLD_SSP245_annualAnom_19852014_20702099_JAS.nc';
readncfile
lon = lon - 360;
[~,~,MLD_anom.jas] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);
fname = 'MLD_SSP245_annualAnom_19852014_20702099_OND.nc';
readncfile
lon = lon - 360;
[~,~,MLD_anom.ond] = clipData(anomaly,lon,lat,lon_min,lon_max,lat_min,lat_max);

MLD_anom.lon_interp = sst.lon; MLD_anom.lat_interp = sst.lat; 
MLD_anom = interpolateData(MLD_anom); % Interpolate

% Remove land points
SST_anom.mask = sst.mask;
SST_anom.land = sst.land;
SST_anom.year_interp(SST_anom.land == 1) = NaN;
SST_anom.jfm_interp(SST_anom.land == 1) = NaN;
SST_anom.amj_interp(SST_anom.land == 1) = NaN;
SST_anom.jas_interp(SST_anom.land == 1) = NaN;
SST_anom.ond_interp(SST_anom.land == 1) = NaN;

MLD_anom.mask = mld.mask;
MLD_anom.land = mld.land;
MLD_anom.year_interp(MLD_anom.land == 1) = NaN;
MLD_anom.jfm_interp(MLD_anom.land == 1) = NaN;
MLD_anom.amj_interp(MLD_anom.land == 1) = NaN;
MLD_anom.jas_interp(MLD_anom.land == 1) = NaN;
MLD_anom.ond_interp(MLD_anom.land == 1) = NaN;

%% Calculate SST and MLD Projections using Anomalies

SST_proj = struct(); SST_proj.lon = sst.lon; SST_proj.lat = sst.lat; SST_proj.mask = sst.mask; SST_proj.land = sst.land;
MLD_proj = struct(); MLD_proj.lon = mld.lon; MLD_proj.lat = mld.lat; MLD_proj.mask = mld.mask; MLD_proj.land = mld.land;

for i = 1:3
    SST_proj.seas(:,:,i) = sst.seas(:,:,i) + SST_anom.jfm_interp;
    MLD_proj.seas(:,:,i) = mld.seas(:,:,i) + MLD_anom.jfm_interp;
end
for i = 4:6
    SST_proj.seas(:,:,i) = sst.seas(:,:,i) + SST_anom.amj_interp;
    MLD_proj.seas(:,:,i) = mld.seas(:,:,i) + MLD_anom.amj_interp;
end
for i = 7:9
    SST_proj.seas(:,:,i) = sst.seas(:,:,i) + SST_anom.jas_interp;
    MLD_proj.seas(:,:,i) = mld.seas(:,:,i) + MLD_anom.jas_interp;
end
for i = 10:12
    SST_proj.seas(:,:,i) = sst.seas(:,:,i) + SST_anom.ond_interp;
    MLD_proj.seas(:,:,i) = mld.seas(:,:,i) + MLD_anom.ond_interp;
end

% Creating new structs for projections to run simulations
PROJECTIONS = struct();

PROJECTIONS.sst = prepStruct(SST_proj,mld,fish,prod,"seas"); % Only SST projection
PROJECTIONS.mld = prepStruct(sst,MLD_proj,fish,prod,"seas"); % Only MLD projection
PROJECTIONS.sstmld = prepStruct(SST_proj,MLD_proj,fish,prod,"seas"); % SST and MLD projection

clear ans cb fname histclim histstddev varratio t lon_min lon_max lat_min lat_max i ll

%% Run simulations under different scenarios

% Run simulations
SIMULATIONS.present = FISH_run_simulation(sst,fish2d,fish,"swim");
SIMULATIONS.sst = FISH_run_simulation(PROJECTIONS.sst,fish2d,fish,"swim");
SIMULATIONS.mld = FISH_run_simulation(PROJECTIONS.mld,fish2d,fish,"swim");
SIMULATIONS.sstmld = FISH_run_simulation(PROJECTIONS.sstmld,fish2d,fish,"swim");

SIMULATIONS.present = calculate_seas(SIMULATIONS.present);
SIMULATIONS.sst = calculate_seas(SIMULATIONS.sst);
SIMULATIONS.mld = calculate_seas(SIMULATIONS.mld);
SIMULATIONS.sstmld = calculate_seas(SIMULATIONS.sstmld);

% Create grid for density maps
lon_edges = linspace(min(sst.lon(:,1)), max(sst.lon(:,1)), 100);
lat_edges = linspace(min(sst.lat(1,:)), max(sst.lat(1,:)), 100);
lon_centers = (lon_edges(2:end) + lon_edges(1:end-1))/2;
lat_centers = (lat_edges(2:end) + lat_edges(1:end-1))/2;

% Calculate density maps
[SIMULATIONS.present.counts, ~, ~] = histcounts2(SIMULATIONS.present.X, SIMULATIONS.present.Y, lon_edges, lat_edges);
[SIMULATIONS.sst.counts, ~, ~] = histcounts2(SIMULATIONS.sst.X, SIMULATIONS.sst.Y, lon_edges, lat_edges);
[SIMULATIONS.mld.counts, ~, ~] = histcounts2(SIMULATIONS.mld.X, SIMULATIONS.mld.Y, lon_edges, lat_edges);
[SIMULATIONS.sstmld.counts, ~, ~] = histcounts2(SIMULATIONS.sstmld.X, SIMULATIONS.sstmld.Y, lon_edges, lat_edges);

% Calculate counts for each scenario
SIMULATIONS.present.counts0 = SIMULATIONS.present.counts'; 
SIMULATIONS.sst.counts0 = SIMULATIONS.sst.counts'; 
SIMULATIONS.mld.counts0 = SIMULATIONS.mld.counts'; 
SIMULATIONS.sstmld.counts0 = SIMULATIONS.sstmld.counts'; 

SIMULATIONS.present.counts(SIMULATIONS.present.counts == 0) = NaN; SIMULATIONS.present.counts = SIMULATIONS.present.counts';
SIMULATIONS.sst.counts(SIMULATIONS.sst.counts == 0) = NaN; SIMULATIONS.sst.counts = SIMULATIONS.sst.counts';
SIMULATIONS.mld.counts(SIMULATIONS.mld.counts == 0) = NaN; SIMULATIONS.mld.counts = SIMULATIONS.mld.counts';
SIMULATIONS.sstmld.counts(SIMULATIONS.sstmld.counts == 0) = NaN; SIMULATIONS.sstmld.counts = SIMULATIONS.sstmld.counts';

% % Calculate difference between projection and present scenarios
% SIMULATIONS.sst.diff_abs = SIMULATIONS.sst.counts - SIMULATIONS.present.counts; SIMULATIONS.sst.diff_abs(SIMULATIONS.sst.diff_abs == 0) = NaN;
% SIMULATIONS.mld.diff_abs = SIMULATIONS.mld.counts - SIMULATIONS.present.counts; SIMULATIONS.mld.diff_abs(SIMULATIONS.mld.diff_abs == 0) = NaN;
% SIMULATIONS.sstmld.diff_abs = SIMULATIONS.sstmld.counts - SIMULATIONS.present.counts; SIMULATIONS.sstmld.diff_abs(SIMULATIONS.sstmld.diff_abs == 0) = NaN;
% 
% SIMULATIONS.sst.diff_ratio = (SIMULATIONS.sst.counts0 - SIMULATIONS.present.counts0)./SIMULATIONS.present.counts0*100;
% SIMULATIONS.mld.diff_ratio = (SIMULATIONS.mld.counts0 - SIMULATIONS.present.counts0)./SIMULATIONS.present.counts0*100;
% SIMULATIONS.sstmld.diff_ratio = (SIMULATIONS.sstmld.counts0 - SIMULATIONS.present.counts0)./SIMULATIONS.present.counts0*100;

% SST scenario difference map
SIMULATIONS.sst.diff_map = SIMULATIONS.sst.counts0 - SIMULATIONS.present.counts0;  % Difference
SIMULATIONS.sst.diff_map(SIMULATIONS.sst.counts0 > 0 & SIMULATIONS.present.counts0 == 0) = 99999;   % Assign 999 where only modelA predicts
SIMULATIONS.sst.diff_map(SIMULATIONS.sst.counts0 == 0 & SIMULATIONS.present.counts0 > 0) = -99999;
SIMULATIONS.sst.diff_map(SIMULATIONS.sst.counts0 == 0 & SIMULATIONS.present.counts0 == 0) = NaN;

% MLD scenario difference map
SIMULATIONS.mld.diff_map = SIMULATIONS.mld.counts0 - SIMULATIONS.present.counts0;  % Difference
SIMULATIONS.mld.diff_map(SIMULATIONS.mld.counts0 > 0 & SIMULATIONS.present.counts0 == 0) = 99999;   % Assign 999 where only modelA predicts
SIMULATIONS.mld.diff_map(SIMULATIONS.mld.counts0 == 0 & SIMULATIONS.present.counts0 > 0) = -99999;
SIMULATIONS.mld.diff_map(SIMULATIONS.mld.counts0 == 0 & SIMULATIONS.present.counts0 == 0) = NaN;

% SSTMLD scenario difference map
SIMULATIONS.sstmld.diff_map = SIMULATIONS.sstmld.counts0 - SIMULATIONS.present.counts0;  % Difference
SIMULATIONS.sstmld.diff_map(SIMULATIONS.sstmld.counts0 > 0 & SIMULATIONS.present.counts0 == 0) = 99999;   % Assign 999 where only modelA predicts
SIMULATIONS.sstmld.diff_map(SIMULATIONS.sstmld.counts0 == 0 & SIMULATIONS.present.counts0 > 0) = -99999;
SIMULATIONS.sstmld.diff_map(SIMULATIONS.sstmld.counts0 == 0 & SIMULATIONS.present.counts0 == 0) = NaN;

sst_sstmld_diff_map = SIMULATIONS.sstmld.counts0 - SIMULATIONS.sst.counts0;  % Difference
sst_sstmld_diff_map(SIMULATIONS.sstmld.counts0 > 0 & SIMULATIONS.sst.counts0 == 0) = 99999;   % Assign 999 where only modelA predicts
sst_sstmld_diff_map(SIMULATIONS.sstmld.counts0 == 0 & SIMULATIONS.sst.counts0 > 0) = -99999;
sst_sstmld_diff_map(SIMULATIONS.sstmld.counts0 == 0 & SIMULATIONS.sst.counts0 == 0) = NaN;

%% Calculate % albacore in CCLME and experiencing MLD <30m

sst_lon_list = sst.lon(:);
sst_lat_list = sst.lat(:);
lme_mask_list = sst.lme_mask(:);
sst_coords = [sst_lon_list, sst_lat_list];

simulations = {SIMULATIONS.present, SIMULATIONS.sst, SIMULATIONS.mld, SIMULATIONS.sstmld};
for i = 1:length(simulations)
    simulation = simulations{i};
    simulation.percent_inside_ccs = zeros(12, 1);
    for imon = 1:12
        in = find(simulation.MON == imon);
        X_positions = simulation.X(in);
        Y_positions = simulation.Y(in);
        
        % Find nearest grid points in the mask for each position
        albacore_coords = [X_positions, Y_positions];
        idx = knnsearch(sst_coords, albacore_coords);
        
        % Use the indices to get mask values for each position
        in_mask = lme_mask_list(idx);
        
        % Calculate percentage inside CCS mask
        simulation.percent_inside_ccs(imon) = 100 * sum(~isnan(in_mask)) / length(in_mask);
    end
    simulations{i} = simulation; % Save the updated struct back
end
SIMULATIONS.present = simulations{1}; SIMULATIONS.sst = simulations{2}; SIMULATIONS.mld = simulations{3}; SIMULATIONS.sstmld = simulations{4};

%% Split simulations by season (density map for reference scenario and difference maps for future scenarios)

% Define the 'imon' matrix with specific month indices

clear imon

imon(:,1) = [1 2 3 0 0 0];
imon(:,2) = [4 5 6 0 0 0];
imon(:,3) = [7 8 9 0 0 0];
imon(:,4) = [10 11 12 0 0 0];

% Indices for each season
jfm = ismember(SIMULATIONS.present.MON, imon(:,1));
amj = ismember(SIMULATIONS.present.MON, imon(:,2));
jas = ismember(SIMULATIONS.present.MON, imon(:,3));
ond = ismember(SIMULATIONS.present.MON, imon(:,4));

% Heatmaps for present scenario
[SIMULATIONS.present.byseason.jfm, ~, ~] = histcounts2(SIMULATIONS.present.X(jfm), SIMULATIONS.present.Y(jfm), lon_edges, lat_edges);
SIMULATIONS.present.byseason.jfm = SIMULATIONS.present.byseason.jfm';
[SIMULATIONS.present.byseason.amj, ~, ~] = histcounts2(SIMULATIONS.present.X(amj), SIMULATIONS.present.Y(amj), lon_edges, lat_edges);
SIMULATIONS.present.byseason.amj = SIMULATIONS.present.byseason.amj'; 
[SIMULATIONS.present.byseason.jas, ~, ~] = histcounts2(SIMULATIONS.present.X(jas), SIMULATIONS.present.Y(jas), lon_edges, lat_edges);
SIMULATIONS.present.byseason.jas = SIMULATIONS.present.byseason.jas'; 
[SIMULATIONS.present.byseason.ond, ~, ~] = histcounts2(SIMULATIONS.present.X(ond), SIMULATIONS.present.Y(ond), lon_edges, lat_edges);
SIMULATIONS.present.byseason.ond = SIMULATIONS.present.byseason.ond';

% Heatmaps for sst scenario
[SIMULATIONS.sst.byseason.jfm, ~, ~] = histcounts2(SIMULATIONS.sst.X(jfm), SIMULATIONS.sst.Y(jfm), lon_edges, lat_edges);
SIMULATIONS.sst.byseason.jfm = SIMULATIONS.sst.byseason.jfm'; 
[SIMULATIONS.sst.byseason.amj, ~, ~] = histcounts2(SIMULATIONS.sst.X(amj), SIMULATIONS.sst.Y(amj), lon_edges, lat_edges);
SIMULATIONS.sst.byseason.amj = SIMULATIONS.sst.byseason.amj'; 
[SIMULATIONS.sst.byseason.jas, ~, ~] = histcounts2(SIMULATIONS.sst.X(jas), SIMULATIONS.sst.Y(jas), lon_edges, lat_edges);
SIMULATIONS.sst.byseason.jas = SIMULATIONS.sst.byseason.jas'; 
[SIMULATIONS.sst.byseason.ond, ~, ~] = histcounts2(SIMULATIONS.sst.X(ond), SIMULATIONS.sst.Y(ond), lon_edges, lat_edges);
SIMULATIONS.sst.byseason.ond = SIMULATIONS.sst.byseason.ond'; 

% Heatmaps for mld scenario
[SIMULATIONS.mld.byseason.jfm, ~, ~] = histcounts2(SIMULATIONS.mld.X(jfm), SIMULATIONS.mld.Y(jfm), lon_edges, lat_edges);
SIMULATIONS.mld.byseason.jfm = SIMULATIONS.mld.byseason.jfm';
[SIMULATIONS.mld.byseason.amj, ~, ~] = histcounts2(SIMULATIONS.mld.X(amj), SIMULATIONS.mld.Y(amj), lon_edges, lat_edges);
SIMULATIONS.mld.byseason.amj = SIMULATIONS.mld.byseason.amj';
[SIMULATIONS.mld.byseason.jas, ~, ~] = histcounts2(SIMULATIONS.mld.X(jas), SIMULATIONS.mld.Y(jas), lon_edges, lat_edges);
SIMULATIONS.mld.byseason.jas = SIMULATIONS.mld.byseason.jas'; 
[SIMULATIONS.mld.byseason.ond, ~, ~] = histcounts2(SIMULATIONS.mld.X(ond), SIMULATIONS.mld.Y(ond), lon_edges, lat_edges);
SIMULATIONS.mld.byseason.ond = SIMULATIONS.mld.byseason.ond';

% Heatmaps for sst+mld scenario
[SIMULATIONS.sstmld.byseason.jfm, ~, ~] = histcounts2(SIMULATIONS.sstmld.X(jfm), SIMULATIONS.sstmld.Y(jfm), lon_edges, lat_edges);
SIMULATIONS.sstmld.byseason.jfm = SIMULATIONS.sstmld.byseason.jfm';
[SIMULATIONS.sstmld.byseason.amj, ~, ~] = histcounts2(SIMULATIONS.sstmld.X(amj), SIMULATIONS.sstmld.Y(amj), lon_edges, lat_edges);
SIMULATIONS.sstmld.byseason.amj = SIMULATIONS.sstmld.byseason.amj';
[SIMULATIONS.sstmld.byseason.jas, ~, ~] = histcounts2(SIMULATIONS.sstmld.X(jas), SIMULATIONS.sstmld.Y(jas), lon_edges, lat_edges);
SIMULATIONS.sstmld.byseason.jas = SIMULATIONS.sstmld.byseason.jas';
[SIMULATIONS.sstmld.byseason.ond, ~, ~] = histcounts2(SIMULATIONS.sstmld.X(ond), SIMULATIONS.sstmld.Y(ond), lon_edges, lat_edges);
SIMULATIONS.sstmld.byseason.ond = SIMULATIONS.sstmld.byseason.ond';

% Calculate difference maps between future scenario and present scenario
SIMULATIONS = season_map_diff(SIMULATIONS,"sst");
SIMULATIONS = season_map_diff(SIMULATIONS,"mld");
SIMULATIONS = season_map_diff(SIMULATIONS,"sstmld");

% Turn 0s to NaN for mapping
SIMULATIONS.present.byseason.jfm(SIMULATIONS.present.byseason.jfm == 0) = NaN; 
SIMULATIONS.present.byseason.amj(SIMULATIONS.present.byseason.amj == 0) = NaN; 
SIMULATIONS.present.byseason.jas(SIMULATIONS.present.byseason.jas == 0) = NaN; 
SIMULATIONS.present.byseason.ond(SIMULATIONS.present.byseason.ond == 0) = NaN; 

SIMULATIONS.sst.byseason.jfm(SIMULATIONS.sst.byseason.jfm == 0) = NaN; 
SIMULATIONS.sst.byseason.amj(SIMULATIONS.sst.byseason.amj == 0) = NaN;
SIMULATIONS.sst.byseason.jas(SIMULATIONS.sst.byseason.jas == 0) = NaN; 
SIMULATIONS.sst.byseason.ond(SIMULATIONS.sst.byseason.ond == 0) = NaN;

SIMULATIONS.mld.byseason.jfm(SIMULATIONS.mld.byseason.jfm == 0) = NaN; 
SIMULATIONS.mld.byseason.amj(SIMULATIONS.mld.byseason.amj == 0) = NaN;
SIMULATIONS.mld.byseason.jas(SIMULATIONS.mld.byseason.jas == 0) = NaN; 
SIMULATIONS.mld.byseason.ond(SIMULATIONS.mld.byseason.ond == 0) = NaN;

SIMULATIONS.sstmld.byseason.jfm(SIMULATIONS.sstmld.byseason.jfm == 0) = NaN; 
SIMULATIONS.sstmld.byseason.amj(SIMULATIONS.sstmld.byseason.amj == 0) = NaN;
SIMULATIONS.sstmld.byseason.jas(SIMULATIONS.sstmld.byseason.jas == 0) = NaN; 
SIMULATIONS.sstmld.byseason.ond(SIMULATIONS.sstmld.byseason.ond == 0) = NaN;

datapoints = [304,485,648,235,255,252,279,664,413,274,641,626];
start_idx = zeros(12,10); start_idx(1,1) = 1;
end_idx = zeros(12,10);
for i = 2:length(datapoints)
    start_idx(i,1) = start_idx(i-1) + datapoints(i-1);
end
for i = 1:length(datapoints)
    end_idx(i,1) = start_idx(i)-1 + datapoints(i);
end
for i = 2:10
    start_idx(:,i) = start_idx(:,i-1) + sum(datapoints);
    end_idx(:,i) = end_idx(:,i-1) + sum(datapoints);
end

for i = 1:12
    first(i) = datetime(fish2d(i).time(1), 'ConvertFrom', 'datenum');
    sim_dates{i} = datetime(fish2d(i).time(1:end-1), 'ConvertFrom', 'datenum');
    obs_dates{i} = datetime(fish2d(i).time, 'ConvertFrom', 'datenum');
end

%% Calculate entry, exit, time in CCS, and migration distance stats for simulated tracks

exits = struct();
entries = struct();
migration_dist = struct();
timeCCS = struct();

scenarios = {"present","sst","mld","sstmld"};
for sims = 1:length(scenarios)
    scen = scenarios{sims};
    for reps = 1:10
        for fish_idx = 1:12

            index_start = start_idx(fish_idx,reps);
            index_end = end_idx(fish_idx,reps);
    
            X_positions = SIMULATIONS.(scen).X(index_start:index_end);
            Y_positions = SIMULATIONS.(scen).Y(index_start:index_end);
    
            % Calculate migration distance
            [~, idx] = min(X_positions); % idx is the index of the most western point
            mig_lon = X_positions(idx);
            mig_lat = Y_positions(idx);
            migration_dist.(scen).mean(fish_idx,reps) = round(haversine(X_positions(1),Y_positions(1),mig_lon,mig_lat)/1000);
    
            isInside = zeros(datapoints(fish_idx),1);
            for i = 1:length(isInside)
                [~, lat_idx] = min(abs(sst.lat(1,:) - Y_positions(i)));  % Index of closest latitude
                [~, lon_idx] = min(abs(sst.lon(:,1) - X_positions(i)));  % Index of closest longitude
                
                % Check if the closest point in the meshgrid is inside the CCS mask
                if sst.lme_mask(lon_idx,lat_idx) == 1
                    isInside(i) = 1;  % Inside the CCS
                else
                    isInside(i) = 0;  % Outside the CCS
                end
            end
            timeCCS.(scen).idx{fish_idx,reps} = isInside;
    
            leave_all = [];
            enter_all = [];
    
            i = 2;
            
            while i <= length(isInside) - 1
                if isInside(i-1) == 1 && isInside(i) == 0
                    leave_all(end+1) = i-1;
                end
                if isInside(i-1) == 0 && isInside(i) == 1
                    enter_all(end+1) = i-1;
                end
                i = i + 1;
            end
    
            min_days = 30;
            true_exit_idx = leave_all;
            true_entry_idx = enter_all;

            marker = 0;
    
            for i = length(leave_all):-1:1 % Iterate in reverse to avoid indexing issues
                exit_idx = leave_all(i);
                if any(enter_all > exit_idx & enter_all <= exit_idx + min_days)
                    marker = 1;
                end
                % if any(enter_all < exit_idx & enter_all >= exit_idx - 10)
                %     marker = 1;
                % end
                if any(leave_all > exit_idx & leave_all <= exit_idx + min_days)
                    marker = 1; % Remove this exit
                end
                if marker == 1
                    true_exit_idx(i) = []; % Remove this exit
                end
                marker = 0;
            end
            
            % Check entries: Remove entries that have an exit within the next 30 days
            for i = length(enter_all):-1:1 % Iterate in reverse to avoid indexing issues
                entry_idx = enter_all(i);
                if any(leave_all < entry_idx & leave_all >= entry_idx - min_days)
                    marker = 1;
                end
                % if any(leave_all > entry_idx & leave_all <= entry_idx + 10)
                %     marker = 1;
                % end
                if any(enter_all < entry_idx & enter_all >= entry_idx - min_days)
                    marker = 1;
                end
                if marker == 1
                    true_entry_idx(i) = []; % Remove this exit
                end
                marker = 0;
            end

            % If any true entries and exits are within 5 days --> remove
            % both
            i = 1;
            while i <= length(true_entry_idx)
                true_entry = true_entry_idx(i);

                % Find the corresponding exit within 10 days
                close_exit_idx = find(true_exit_idx > true_entry & true_exit_idx <= true_entry + 5, 1);

                if ~isempty(close_exit_idx)
                    % Remove both the entry and exit indices
                    true_exit_idx(close_exit_idx) = [];
                    true_entry_idx(i) = [];

                    % Restart the loop to ensure we correctly process remaining values
                    i = 1;
                else
                    i = i + 1;
                end
            end
    
            exits.(scen).idx{fish_idx,reps} = true_exit_idx;
            entries.(scen).idx{fish_idx,reps} = true_entry_idx;

            dates = fish2d(fish_idx).time(1:end-1);
            day_of_year = datenum(datestr(dates, 'yyyy-mm-dd')) - datenum(datestr(dates, 'yyyy-01-01'));
            exits.(scen).date{fish_idx,reps} = day_of_year(true_exit_idx);
            entries.(scen).date{fish_idx,reps} = day_of_year(true_entry_idx);

            % Calculate mean entry and exit for fish that have multiple
            % using radians
            if ~isempty(exits.(scen).date{fish_idx,reps})
                rad = exits.(scen).date{fish_idx,reps}/365*2*pi;
                x = cos(rad);
                y = sin(rad);
                angle = atan2(mean(y), mean(x));
                exits.(scen).mean(fish_idx,reps) = mod(angle,2*pi)/(2*pi)*365;
            else
                exits.(scen).mean(fish_idx,reps) = NaN;
            end
            if ~isempty(entries.(scen).date{fish_idx,reps})
                rad = entries.(scen).date{fish_idx,reps}/365*2*pi;
                x = cos(rad);
                y = sin(rad);
                angle = atan2(mean(y), mean(x));
                entries.(scen).mean(fish_idx,reps) = mod(angle,2*pi)/(2*pi)*365;
            else
                entries.(scen).mean(fish_idx,reps) = NaN;
            end
            
            if length(dates) >= 365
                timeCCS.(scen).mean(fish_idx,reps) = sum(isInside(1:365));
            else 
                timeCCS.(scen).mean(fish_idx,reps) = NaN;
            end

            if ~isempty(true_exit_idx) && ~isempty(true_entry_idx)
                timeCCS.(scen).mean2(fish_idx,reps) = mod(exits.(scen).mean(fish_idx,reps) - entries.(scen).mean(fish_idx,reps), 365);
            else
                timeCCS.(scen).mean2(fish_idx,reps) = NaN;
            end
        end
    end
end

%% Calculate entry, exit, time in CCS, and migration distance stats for tags

for fish_idx = 1:12
    index_start = start_idx(fish_idx,reps);
    index_end = end_idx(fish_idx,reps);
    
    X_positions = fish2d(fish_idx).lon;
    Y_positions = fish2d(fish_idx).lat;
    
    % Calculate migration distance
    [~, idx] = min(X_positions); % idx is the index of the most western point
    mig_lon = X_positions(idx);
    mig_lat = Y_positions(idx);
    migration_dist.stats.mean(fish_idx) = round(haversine(X_positions(1),Y_positions(1),mig_lon,mig_lat)/1000);
    
    isInside = zeros(datapoints(fish_idx),1);
    for i = 1:length(isInside)
        [~, lat_idx] = min(abs(sst.lat(1,:) - Y_positions(i)));  % Index of closest latitude
        [~, lon_idx] = min(abs(sst.lon(:,1) - X_positions(i)));  % Index of closest longitude
        
        % Check if the closest point in the meshgrid is inside the CCS mask
        if sst.lme_mask(lon_idx,lat_idx) == 1
            isInside(i) = 1;  % Inside the CCS
        else
            isInside(i) = 0;  % Outside the CCS
        end
    end
    timeCCS.stats.idx{fish_idx} = isInside;
    
    leave_all = [];
    enter_all = [];
    
    i = 2;
    
    while i <= length(isInside) - 1
        if isInside(i-1) == 1 && isInside(i) == 0
            leave_all(end+1) = i-1;
        end
        if isInside(i-1) == 0 && isInside(i) == 1
            enter_all(end+1) = i-1;
        end
        i = i + 1;
    end
    
    min_days = 30;
    true_exit_idx = leave_all;
    true_entry_idx = enter_all;
    
    for i = length(leave_all):-1:1 % Iterate in reverse to avoid indexing issues
        exit_idx = leave_all(i);
        if any(enter_all > exit_idx & enter_all <= exit_idx + min_days)
            marker = 1;
        end
        % if any(enter_all < exit_idx & enter_all >= exit_idx - 10)
        %     marker = 1;
        % end
        if any(leave_all > exit_idx & leave_all <= exit_idx + min_days)
            marker = 1; % Remove this exit
        end
        if marker == 1
            true_exit_idx(i) = []; % Remove this exit
        end
        marker = 0;
    end
    
    % Check entries: Remove entries that have an exit within the next 15 days
    for i = length(enter_all):-1:1 % Iterate in reverse to avoid indexing issues
        entry_idx = enter_all(i);
        if any(leave_all < entry_idx & leave_all >= entry_idx - min_days)
            marker = 1;
        end
        % if any(leave_all > entry_idx & leave_all <= entry_idx + 10)
        %     marker = 1;
        % end
        if any(enter_all < entry_idx & enter_all >= entry_idx - min_days)
            marker = 1;
        end
        if marker == 1
            true_entry_idx(i) = []; % Remove this exit
        end
        marker = 0;
    end
    
    exits.stats.idx{fish_idx} = true_exit_idx;
    entries.stats.idx{fish_idx} = true_entry_idx;
    
    dates = fish2d(fish_idx).time;
    day_of_year = datenum(datestr(dates, 'yyyy-mm-dd')) - datenum(datestr(dates, 'yyyy-01-01'));
    exits.stats.date{fish_idx} = day_of_year(true_exit_idx);
    entries.stats.date{fish_idx} = day_of_year(true_entry_idx);

    if ~isempty(exits.stats.date{fish_idx})
        rad = exits.stats.date{fish_idx}/365*2*pi;
        x = cos(rad);
        y = sin(rad);
        angle = atan2(mean(y), mean(x));
        exits.stats.mean(fish_idx) = mod(angle,2*pi)/(2*pi)*365;
    end
    if ~isempty(entries.stats.date{fish_idx})
        rad = entries.stats.date{fish_idx}/365*2*pi;
        x = cos(rad);
        y = sin(rad);
        angle = atan2(mean(y), mean(x));
        entries.stats.mean(fish_idx) = mod(angle,2*pi)/(2*pi)*365;
    end
            
    if length(dates) >= 365
        timeCCS.stats.mean(fish_idx) = sum(isInside(1:365));
    else 
        timeCCS.stats.mean(fish_idx) = NaN;
    end

    if ~isempty(true_exit_idx) && ~isempty(true_entry_idx)
        timeCCS.stats.mean2(fish_idx) = mod(exits.stats.mean(fish_idx) - entries.stats.mean(fish_idx), 365);
    else
        timeCCS.stats.mean2(fish_idx) = NaN;
    end 
end

%% Process entry, exit, time in CCS, and migration distance data

scenarios = {"present","sst","mld","sstmld"};
migration_dist_data = nan(120, length(scenarios));
timeCCS_data = nan(120, length(scenarios));
timeCCS_data2 = nan(120, length(scenarios));

entry_data = nan(120, length(scenarios));
entry_data_mean = nan(length(scenarios),1);
exit_data = nan(120, length(scenarios));
exit_data_mean = nan(length(scenarios),1);

for sims = 1:length(scenarios)
    scen = scenarios{sims};

    % Non circular data
    migration_dist_data(:,sims) = migration_dist.(scen).mean(:);
    timeCCS_data(:,sims) = timeCCS.(scen).mean(:);
    timeCCS_data2(:,sims) = timeCCS.(scen).mean2(:);

    % Circular data
    [exit_data(:,sims), exit_data_mean(sims)] = circular_doy(exits.(scen).mean(:));
    [entry_data(:,sims), entry_data_mean(sims)] = circular_doy(entries.(scen).mean(:));
end

scenarios = {"present","sst","mld","sstmld"};
timeCCS.processed.all = nan(120, length(scenarios));
timeCCS.processed.all2 = nan(120, length(scenarios));
migration_dist.processed.all = nan(120, length(scenarios));
entries.processed.all = nan(120, length(scenarios));
exits.processed.all = nan(120, length(scenarios));

for sims = 1:length(scenarios)
    scen = scenarios{sims};

    migration_dist.processed.all(:,sims) = migration_dist.(scen).mean(:);
    timeCCS.processed.all(:,sims) = timeCCS.(scen).mean(:);
    timeCCS.processed.all2(:,sims) = timeCCS.(scen).mean2(:);
    entries.processed.all(:,sims) = entries.(scen).mean(:);
    exits.processed.all(:,sims) = exits.(scen).mean(:);
end

[~,entries.processed.mean] = circular_doy(entries.processed.all(:));
[~,exits.processed.mean] = circular_doy(exits.processed.all(:));

% Shift datasets to global reference
for i = 1:length(scenarios)
    entries.processed.shifted(:,i) = mod(entries.processed.all(:,i) - entries.processed.mean + 182.5, 365);
    exits.processed.shifted(:,i) = mod(exits.processed.all(:,i) - exits.processed.mean + 182.5, 365);
end

%% Calculate change in center of gravity between reference and future scenario

R = 6371;

% Preallocate
cog_shift_km = zeros(12,1);
cog_bearing_deg = zeros(12,1);

lat_present = SIMULATIONS.present.latseas;
lon_present = SIMULATIONS.present.lonseas;
lat_future = SIMULATIONS.sstmld.latseas;
lon_future = SIMULATIONS.sstmld.lonseas;

for m = 1:12
    % Convert degrees to radians
    lat1 = deg2rad(lat_present(m));
    lon1 = deg2rad(lon_present(m));
    lat2 = deg2rad(lat_future(m));
    lon2 = deg2rad(lon_future(m));
    
    % Haversine distance
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;
    
    a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    cog_shift_km(m) = R * c;
    
    % Bearing calculation
    y = sin(dlon) * cos(lat2);
    x = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dlon);
    theta = atan2(y, x);  % in radians
    bearing_deg = mod(rad2deg(theta), 360);  % convert to 0–360°
    
    cog_bearing_deg(m) = bearing_deg;
end

SIMULATIONS.sstmld.cog_shift = cog_shift_km;
SIMULATIONS.sstmld.cog_shift_bearing = cog_bearing_deg;

%% Save all structs

save('processed_data/mld.mat', 'mld')
save('processed_data/sst.mat', 'sst','-v7.3')
save('processed_data/MLD_anom.mat', 'MLD_anom')
save('processed_data/SST_anom.mat', 'SST_anom')
save('processed_data/MLD_proj.mat', 'MLD_proj')
save('processed_data/SST_proj.mat', 'SST_proj')
save('processed_data/PROJECTIONS.mat', 'PROJECTIONS','-v7.3')
save('processed_data/fish.mat', 'fish')
save('processed_data/fish2d.mat', 'fish2d')
save('processed_data/fish_vertical.mat', 'fish_vertical')
save('processed_data/mig_fish.mat', 'mig_fish')
save('processed_data/indiv_fish.mat', 'indiv_fish')
save('processed_data/SIMULATIONS.mat', 'SIMULATIONS')
save('processed_data/lat_centers.mat', 'lat_centers')
save('processed_data/lon_centers.mat', 'lon_centers')
save('processed_data/exits.mat', 'exits')
save('processed_data/entries.mat', 'entries')
save('processed_data/migration_dist.mat', 'migration_dist')
save('processed_data/timeCCS.mat', 'timeCCS')

