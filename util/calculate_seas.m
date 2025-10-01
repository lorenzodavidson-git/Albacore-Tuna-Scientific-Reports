function updated_struct = calculate_seas(struct)
    for imon = 1:12
        in = find(struct.MON == imon);
        struct.useas(imon) = nanmean(struct.U(in));
        struct.vseas(imon) = nanmean(struct.V(in));
        struct.lonseas(imon) = nanmean(struct.X(in));
        struct.latseas(imon) = nanmean(struct.Y(in));
        struct.mldseas(imon) = nanmean(struct.MLD(in));
    end
    updated_struct = struct;