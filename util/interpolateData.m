function [struct] = interpolateData(struct)

    fields = fieldnames(struct);
    for i = 3:7
        data = struct.(fields{i});
        newfield = [fields{i}, '_interp'];
        interpolant = griddedInterpolant(struct.lon, struct.lat, data, 'linear');
        struct.(newfield) = interpolant(struct.lon_interp,  struct.lat_interp);
    end
end