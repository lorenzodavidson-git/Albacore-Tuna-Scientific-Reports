
function pac_coast(varargin)

if nargin == 0
   color='k';
   varargin{1} = color;
end

load(which('world_coast_data.mat'));
plot(lon,lat,varargin{:})
