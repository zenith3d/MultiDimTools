function out = integrate(obj, dim_name)
    % INTEGRATE Détermine l'intégrale discrète de la variable selon
    % la dimension désirée.
    
    % Détermination de la dimension concernée
    [~, i_objdim, ~] = ...
        intersect(obj.dim_names, dim_name, 'stable');
    
    % Si aucune dimension n'est concérnée ou que la dimension 
    % concernée ne dispose que d'un seul point (singleton),
    % l'intégrale est considérée comme nulle
    if isempty(i_objdim) || (obj.shape(i_objdim) < 2)
        out = obj;
        out.values(:) = 0;
        return;
    end
    
    % Récupération des points correspondants à la dimension 
    % concernée
    pts       = obj.dim_points{i_objdim};
    end_index = obj.shape(i_objdim);
    
    % Initialisation de la table de sortie 
    intC = zeros(obj.shape);
    
    % Définition de la structure 'subs' par défaut
    default_struct.type = '()';
    default_struct.subs = repmat({':'}, 1, obj.n_dims);

    % Intégration discrète de la table (la première valeur de la
    % table est supposée nulle)
    for i = 2:end_index
        SD = default_struct; SD.subs{i_objdim} = i;
        SL = default_struct; SL.subs{i_objdim} = i-1;
        intval = squeeze(subsref(intC,SL)) + ...
                (squeeze(subsref(obj.values,SL)) + ...
                 squeeze(subsref(obj.values,SD))) * ...
                (pts(i) - pts(i-1)) * 0.5;
        intC = subsasgn(intC, SD, intval);
    end
    
    % Création de l'instance de sortie
    out = MultiDimVar(intC, obj.dim_names, obj.dim_points);
    
end