function [clipped_lon, clipped_lat,clipped_matrix] = clipData(matrix, lon, lat, lon_min, lon_max, lat_min, lat_max)
  
    % Find indices within the longitude bounds
    lon_idx = lon >= lon_min & lon <= lon_max;
    
    % Find indices within the latitude bounds
    lat_idx = lat >= lat_min & lat <= lat_max;
    
    % Clip the matrix to the specified bounds
    if ndims(matrix) > 2
        dims = size(matrix);
        for i = 1:dims(3)
            clipped_matrix(:,:,i) = matrix(lon_idx, lat_idx,i);
        end
    else
        clipped_matrix = matrix(lon_idx, lat_idx);
    end
    
    
    % Update latitude and longitude vectors
    lon = lon(lon_idx);
    lat = lat(lat_idx);
    [clipped_lat,clipped_lon] = meshgrid(lat,lon);
end