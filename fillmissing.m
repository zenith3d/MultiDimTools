function out = fillmissing(obj, dim_name, method)
    % FILLMISSING Remplie les entrées manquantes dans le tableau
    % selon la méthode utilisée.
    
    % Détermination de la dimension concernée
    [~, i_objdim, ~] = ...
        intersect(obj.dim_names, dim_name, 'stable');
    
    % Vérification
    if isempty(i_objdim)
        out = obj;
        return;
    end
    
    % Remplissage des valeurs manquantes selon la méthode choisie
    out = obj;
    out.values = fillmissing(obj.values, method, i_objdim, ...
        'EndValues', 'extrap', ...
        'SamplePoints', obj.dim_points{i_objdim});
    
end