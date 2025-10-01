function fish = load_albacore_fish_data


if exist('albacore_fish.mat')
    disp('Loading albacore_fish.mat ... ')
    load albacore_fish.mat
    disp ('Done.')

else

    disp('Processing  NOAA-Albacore-Tag-Data.csv ... ')
    T = readtable('NOAA-Albacore-Tag-Data.csv');
    % initialize fish data structure
    clear fish
    fish.tags = T.tag; %
    fish.year = T.yr;
    fish.month=T.mo;
    fish.day=T.day;
    fish.lon=T.lon360-360;
    fish.lat=T.lat;
    fish.tb=T.tb;
    fish.ta=T.ta;
    fish.light=T.light;
    fish.depth=T.z;
    fish.datenum=datenum(T.dateTime);
    fish.dn=T.dn;
    disp ('Saving the output in albacore_fish.mat')
    save -v7 ../data/albacore_fish.mat fish
    disp ('Done.')
end

%find all the unique tags
fish.tagsid = unique(fish.tags);
fish.numfish = length(fish.tagsid);
