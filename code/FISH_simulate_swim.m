function [x, y, tf, u, v, month] = FISH_simulate_swim(lon, lat, u, v, dt, time, sst, fish)
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
%   mld  - Structure containing mixed layer depth data (not used in this function)
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

% Simulate fish swimming behavior
for t = 1:length(dt)
    
    % Get the current month
    imon = month(t);
    
    % Update velocities
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
    
    % Check if the fish is swimming in acceptable conditions
    % flag = 0; % do not use temperature constraint
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
    
    % Store the temperature experienced by the fish
    tf(t + 1) = temp;
    
    % Generate random velocity for the next time step
    u(t + 1) = a * u(t) + (1 - bu) * u_rand(t);
    v(t + 1) = a * v(t) + (1 - bv) * v_rand(t);
    
    %     % Seasonal rules #1 (works but is prescibing u)
    %     %if ( (ccs_prod(imon) > 0 || food.mld < 30) && food.prod == 0)
    %     if ( (ccs_prod(imon) > 0 || food.mld < 30)  && food.prod == 0)
    %         useas(imon) = fish.u.seas(imon);
    %         fac_u = 0.8;
    %         modified1(t)=1;
    %     end
    %     %if ( (ccs_prod(imon) < 0 || food.mld > 30) && food.prod ~= 0)
    %     if ( ccs_prod(imon) < 0 && food.prod ~= 0)
    % %         useas(imon) = fish.u.seas(imon);
    % %         fac_u = 0.8;
    %          modified2(t)=1;
    %     end
    
    % Seasonal rules #2 (works but is prescibing u)
    
    %     if ( (ccs_prod(imon) > 0 || food.mld < 30) && food.prod == 0)
    %         [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1), food.theta, weight);
    %         u(t+1) = nudgedU;
    %         v(t+1) = nudgedV;
    %         modified1(t)=1;
    %     end
    %     %if ( (ccs_prod(imon) < 0 || food.mld > 30) && food.prod ~= 0)
    
    
    % nudgedU = u(t+1);
    % nudgedV = v(t+1);
    % 
    % 
    % in spring if MLD is too shallow start going back to CCS
    if imon >= 1 && imon <=5 && food.mld < 30 && food.prod == 0
        weight = 0.3;
        [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1),food.theta, weight);
        modified2(t)=1;
    end

    % during late spring and summer when CCS is productive and MLDs are
    % shawllow go to CCS unless you are there already
    if imon > 5 && imon <= 9 && food.prod == 0 && food.mld < 30
        weight = 0.65;
        [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1),food.theta, weight);
        modified2(t)=1;
    end

    % in october as MLD begins to deep offshore search for deeper MLDs > 50
    if imon >= 10 && food.mld < 50
        weight = 0.45;
        [nudgedU, nudgedV] = FISH_NudgeDirection(u(t+1), v(t+1),180, weight);
        modified2(t)=1;
    end


    u(t+1) = nudgedU;
    v(t+1) = nudgedV;
    
end

% Make the output column vectors
u = u(1:end-1)';
v = v(1:end-1)';
x = x';
y = y';
save modified.mat modified1 modified2

