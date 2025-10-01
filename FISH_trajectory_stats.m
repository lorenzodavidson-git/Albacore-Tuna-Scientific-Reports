function f = FISH_trajectory_stats(fish2d,varname, flag)

% analysis of temperature distribution
f.data=[]; f.month=[];

for i = 1:12
    fish2d(i).ta = fish2d(i).ta';
    month = str2num(datestr(fish2d(i).time, 'mm'));
    if flag == 1
        month = month(1:end-1);
    end
    str=['f.data = [f.data; fish2d(i).',varname,'];'];
    eval(str);
    f.month = [f.month; month];
end

for imon = 1:12
    in = find(f.month == imon);
    f.seas(imon) = meanNaN(f.data(in),1);
    f.seas_std(imon) = stdNaN(f.data(in),1);
end

