function [first_entry_idx, first_exit_idx] = find_first_entry_exit(datenum_vector, isInside)

    years = year(datetime(datenum_vector, 'ConvertFrom', 'datenum')); % Extract years
    unique_years = unique(years); % Get unique years

    first_entry_idx = NaN(length(unique_years), 1);
    first_exit_idx = NaN(length(unique_years), 1);

    for y = 1:length(unique_years)
        year_idx = find(years == unique_years(y)); % Indices of this year
        inside_year = isInside(year_idx); % Subset of inside/outside for the year
        
        % Find the first entry (0 → 1 transition)
        for i = 2:length(inside_year)
            if inside_year(i-1) == 0 && inside_year(i) == 1
                first_entry_idx(y) = year_idx(i);
                break;
            end
        end
        
        % Find the first exit (1 → 0 transition after entry)
        if ~isnan(first_entry_idx(y))
            for i = find(year_idx == first_entry_idx(y)):length(inside_year)
                if inside_year(i-1) == 1 && inside_year(i) == 0
                    first_exit_idx(y) = year_idx(i);
                    break;
                end
            end
        end
    end
end