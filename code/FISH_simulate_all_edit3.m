function [x, y, tf, u, v, month, mld_fish] = FISH_simulate_all_edit3(lon, lat, u, v, dt, time, sst, fish,type)
% FISH_simulate_swim: Simulates fish swimming behavior based on velocity, temperature, and constraints
%
% This function simulates the swimming behavior of a fish by updating its position
% and velocity based on given parameters, checking for temperature constraints,
% and ensuring it stays within acceptable conditions.
%
% Inputs:
%   lon  - Initial longitude of the fish
%   lat  - Initial latitude of the fish
%   u    - Initial eastward velocity component
%   v    - Initial northward velocity component
%   dt   - Time step for the simulation
%   time - Time array for the simulation
%   sst  - Structure containing sea surface temperature data and related fields
%   mld  - Structure containing mixed layer depth data
%   
%   type - 1 of 3 options: null (null hypothesis), temp (only temperature
%   constraints), swim (temp constraints + MLD/prod drivers)
%
% Outputs:
%   x  - Updated longitude positions of the fish over time
%   y  - Updated latitude positions of the fish over time
%   tf - Temperature experienced by the fish at each time step
%   u  - Updated eastward velocity component over time
%   v  - Updated northward velocity component over time

% Extract the month from the time array
month = str2num(datestr(time, 'mm'));

% Initialize position arrays
x = lon;
y = lat;

% Velocity parameters
bu = 0.17;
bv = 0.39;
std_u = 0.7;
std_v = 0.3;
a = 0.6097;

% Seasonal cycle of CCS LME
ccs_prod =[ -1.2204   -1.2204 -1.1205 -0.5207 -0.0208 0.5790 1.1788 ...
    1.3787 1.2787 0.6789 -0.1208 -0.8705];

% Initialize velocity arrays
u(:) = 0;
v(:) = 0;
u(1) = randn * std_u;
v(1) = randn * std_v;

% Add prescribed seasonal velocities
useas = fish.u.seas;
vseas = fish.v.seas;
fac_u = 0.8;
fac_v = 0.8;

% remove seasonal
useas(:) = 0;
vseas(:) = 0;
fac_u = 1;
fac_v = 1;
modified1 = u*0;
modified2 = u*0;

% Generate random velocity components
N = length(u);
u_rand = randn(N + 1, 1) * std_u;
v_rand = randn(N + 1, 1) * std_v;

% Initialize temperature experienced by the fish
tf = zeros(size(dt));

% initialize MLD memory array
mld_fish = zeros(size(dt));

% Simulate fish swimming behavior
for t = 1:length(dt)
    
    % Get the current month
    imon = month(t);
    
    % Update velocitites
    u(t) = u(t)*fac_u + useas(imon);
    v(t) = v(t)*fac_v + vseas(imon);
    
    % Update fish position based on current velocity
    x(t + 1) = x(t) + u(t) * dt(t);
    y(t + 1) = y(t) + v(t) * dt(t);
    
    % reset the factors & seasonal velocities
    fac_u = 1;
    fac_v = 1;
    useas(imon) = 0;
    vseas(imon) = 0;

    if type == "null"
        % Check if fish is swimming onto land and reroute
        flag = 0; % use land buffer
        [temp, myflag, theta, p, food] = FISH_position_check(x(t + 1), y(t + 1), ...
        imon, sst, flag);

        if myflag == false
            amp = sqrt(u(t)^2 + v(t)^2);
            u(t) = amp * cos(theta) * sign(p);
            v(t) = amp * sin(theta) * sign(p);
            x(t + 1) = x(t) + u(t) * dt(t);
            y(t + 1) = y(t) + v(t) * dt(t);
        end
    
        % Generate random velocity for the next time step
        u(t + 1) = a * u(t) + (1 - bu) * u_rand(t);
        v(t + 1) = a * v(t) + (1 - bv) * v_rand(t);
    else
        %  Check if the fish is swimming in acceptable conditions
        flag = 1; % use temperature constraint
        [temp, myflag, theta, p, food] = FISH_position_check(x(t + 1), y(t + 1), ...
        imon, sst, flag);
        
        % If conditions are not acceptable, adjust velocity direction
        if myflag == false
            amp = sqrt(u(t)^2 + v(t)^2);
            u(t) = amp * cos(theta) * sign(p);
            v(t) = amp * sin(theta) * sign(p);
            x(t + 1) = x(t) + u(t) * dt(t);
            y(t + 1) = y(t) + v(t) * dt(t);
        end

        % Generate random velocity for the next time step
        u(t + 1) = a * u(t) + (1 - bu) * u_rand(t);
        v(t + 1) = a * v(t) + (1 - bv) * v_rand(t);
    end

    if type == "swim"
        nudgedU = u(t+1);
        nudgedV = v(t+1);
        
        [~,~,~,~,food] = FISH_position_check(x(t + 1), y(t + 1),imon, sst, flag);
        mld_fish(t)= food.mld;

        if nnz(mld_fish) > 14
            prev_mean = mean(mld_fish(t-14:t-8)); % Mean MLD 2 weeks prior
            new_mean = mean(mld_fish(t-7:t)); % Mean MLD in prior week

            % Rule 1: if MLD is deeper than 30m and is deepening, nudge albacore away from CCS
            if food.mld > 30 && food.mld < 70 && prev_mean - 0.25 < new_mean
                weight = 0.4;
                [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1),180, weight);
                modified2(t)=1;
            end

            % Rule 2: If albacore is not in CCS and MLD is shallower than 30m + is shallowing, nudge albacore towards CCS
            if food.mld < 30 && food.prod == 0 && prev_mean >= new_mean - 0.25
                weight = 0.00225*(food.mld-30)^2+0.1;
                weight = min(max(weight, 0), 0.65); % Quadratic fit with MLD 30m = 0 and MLD ~13m or shallower = 0.65
                [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1),food.theta, weight);
                modified2(t)=1;
            end
        else
            if food.mld > 30 && food.mld < 70
                weight = 0.4;
                [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1),180, weight);
                modified2(t)=1;
            end
            if food.mld < 30 && food.prod == 0
                weight = 0.00225*(food.mld-30)^2+0.1;
                weight = min(max(weight, 0), 0.65); % Quadratic fit with MLD 30m = 0 and MLD ~13m or shallower = 0.65
                [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1),food.theta, weight);
                modified2(t)=1;
            end
        end

        u(t+1) = nudgedU;
        v(t+1) = nudgedV;
    end

    % Final storing of temperature and mld experienced by fish
    [temp,~,~,~,food] = FISH_position_check(x(t + 1), y(t + 1),imon, sst, flag);
    tf(t + 1) = temp;
    mld_fish(t)= food.mld;
end
   
% Make the output column vectors
u = u(1:end-1)';
v = v(1:end-1)';
x = x';
y = y';
% save modified.mat modified1 modified2