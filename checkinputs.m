function err = checkinputs(values, dim_names, dim_points)

    % Vérification de la cohérence des entrées
    values_ndims = ndims(values);
    points_ndims = length(dim_names);
    names_ndims  = length(dim_names);
    values_shape = size(values);
    points_shape = cellfun(@length, dim_points);
    
    if values_ndims ~= points_ndims || values_ndims ~= names_ndims
        err = 1;
    elseif values_shape ~= points_shape
        err = 2;
    else
        err = 0x00;
    end
    
    end