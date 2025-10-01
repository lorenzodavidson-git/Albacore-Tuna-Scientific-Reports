
% The GLORYS12V1 product is the CMEMS global ocean eddy-resolving (1/12°
% horizontal resolution, 50 vertical levels) reanalysis covering the
% altimetry (1993 onward).
% https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description

% DATASET ID: cmems_mod_glo_phy_my_0.083deg_P1M-m (1993 - 2021)
% DATASET ID: cmems_mod_glo_phy_myint_0.083deg_P1M-m (2021 - 2023)

% Process the MLD = Mixed Layer Depth from the
% 01/1993 - 01/2021

files = {'cmems_mod_glo_phy_my_0.083deg_P1M-m_1713216547706.nc' ...
    'cmems_mod_glo_phy_my_0.083deg_P1M-m_1713216805414.nc' ...
    'cmems_mod_glo_phy_my_0.083deg_P1M-m_1713217144341.nc' ...
    'cmems_mod_glo_phy_my_0.083deg_P1M-m_1713217618026.nc' ...
    'cmems_mod_glo_phy_my_0.083deg_P1M-m_1713222368918.nc' ...
    'cmems_mod_glo_phy_my_0.083deg_P1M-m_1713222615084.nc'};

time1 = [];
for i=1:length(files)
        file = files{i};
    nc = netcdf(file);
    time1 = [time1 ; nc{'time'}(:)];
end


for i=1:length(files)
    file = files{i};
    ncdisp(file);
    pause
end


clear time year month
n = 0;
for yr = 1999:2021
    for imon = 1:12
        if yr== 2021 & imon == 1
            return
        end
        
        n = n+1;
        time(n) = datenum(yr, imon, 15);
        year(n) = yr;
        month(n) = imon;
    end
end

MLD=[];

for i=2:length(files)
    i
    file = files{i};
    
    nc = netcdf(file);
    lon = nc{'longitude'}(:);
    lat = nc{'latitude'}(:);
    mld = nc{'mlotst'}(:,:,:);
    close(nc);
    
    [lat, lon] = meshgrid(lat, lon);
    %mld = permute(mld, [3 2 1]);
    MLD = [MLD; mld];
end

mld = permute(MLD, [3 2 1]);
figure;
clf
pcolor(lon, lat, mld(:,:,1)); colorbar; shading flat
caxis ([0 150])
set(gca, 'FontSize',16);

for imon = 1:12
    in = find(month == imon);
    mld_seas(:,:,imon) = mean(mld(:,:,in),3);
end


figure;

for imon =1:12
    clf
    pcolor(lon, lat, mld_seas(:,:,imon)); colorbar; shading flat
    caxis ([0 150])
    set(gca, 'FontSize',16);
    title (imon)
    pause
end

clear mld;
mld.lon = lon;
mld.lat = lat;
mld.seas= mld_seas;
mask = mld.seas(:,:,1);
mask(~isnan(mask))=1;
mld.mask = mask;





