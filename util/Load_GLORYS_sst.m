function sst = Load_GLORYS_sst

% The GLORYS12V1 product is the CMEMS global ocean eddy-resolving (1/12° horizontal resolution, 50 vertical levels) reanalysis covering the altimetry (1993 onward).
if ~exist('GLORYS_sst.mat')
    disp('Processing Copernicus SST reanalysis')
% load the ocean reanalysis 
%open link to temperature data
file='https://edilorenzo:!Kimberly31332@my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1M-m'
nc = netcdf(file);
time = nc{'time'} (:)/(24) + datenum('1950-01-01');
lon = nc{'longitude'}(:);
lat = nc{'latitude'}(:);
[lat,lon] = meshgrid(lat, lon);
iswap = find (lon(:,1)>0);
lon2 = [lon(iswap,:)-360 ; lon(1:iswap(1)-1,:)];
depth = nc{'depth'}(:);

sst.time = time;
sst.year = str2num(datestr(sst.time, 'yyyy'));
sst.month = str2num(datestr(sst.time, 'mm'));

I = find (lon2(:,1) > -250 & lon2(:,1) < -100);
J = find (lat(1,:) > 0 & lat(1,:)  < 62);
T_range = find (sst.year >= 2003 & sst.year <= 2013);

sst.lon = lon2(I,J);
sst.lat = lat(I,J);
sst.depth = depth;

% swap the data matrix (need to re-center coordinates over the pacific)
scale = nc{'thetao'}.scale_factor (:);
fillvalue = nc{'thetao'}.FillValue_ (:);
offset = nc{'thetao'}.add_offset(:);



k=0;
for t = T_range'

    % rules for loading temperature
    tic
    temperature = nc{'thetao'}(t,1,:,:);
    in = find (temperature == fillvalue);
    temperature(in) = nan;
    temperature = temperature*scale + offset;
    temperature = temperature';
    temperature = [temperature(iswap,:) ; temperature(1:iswap(1)-1,:)];
    k = k+1;
    sst.data(:,:,k) = temperature(I,J);
    toc
end

sst.time = sst.time(T_range);
sst.year = str2num(datestr(sst.time, 'yyyy'));
sst.month = str2num(datestr(sst.time, 'mm'));
disp ('Saving GLORYS_sst.mat ...')
save -v7 GLORYS_sst.mat sst
disp('Done')
else 
    disp('Loading GLORYS_sst.mat ... ')
    load GLORYS_sst.mat sst
    disp('Done')
end
