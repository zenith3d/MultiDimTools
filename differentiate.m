function out = differentiate(obj, dim_name)
    % DIFFERENTIATE Calcule la dérivée de la variable en fonction 
    % de la dimension désirée.
    
    % Détermination de la dimension concernée
    [~, i_objdim, ~] = ...
        intersect(obj.dim_names, dim_name, 'stable');
    
    % Si aucune dimension n'est concérnée ou que la dimension 
    % concernée ne dispose que d'un seul point (singleton), la 
    % dérivée est nulle
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
    dCdD = zeros(obj.shape);
    
    % Définition de la structure 'subs' par défaut
    default_struct.type = '()';
    default_struct.subs = repmat({':'}, 1, obj.n_dims);

    % Dérivée en debut de table (ordre 1)
    SD = default_struct; SD.subs{i_objdim} = 1;
    SH = default_struct; SH.subs{i_objdim} = 2;
    der = (squeeze(subsref(obj.values,SH)) - ...
           squeeze(subsref(obj.values,SD))) / ...
          (pts(2) - pts(1));
    dCdD = subsasgn(dCdD, SD, der);

    % Dérivée en fin de table (ordre 1)
    SD = default_struct; SD.subs{i_objdim} = end_index;
    SL = default_struct; SL.subs{i_objdim} = end_index-1;
    der = (squeeze(subsref(obj.values,SD)) - ...
           squeeze(subsref(obj.values,SL))) / ...
          (pts(end_index) - pts(end_index-1));
    dCdD = subsasgn(dCdD, SD, der);

    % Dérivée sur le reste de la table (ordre 2 - différence 
    % centrale pondérée)
    for i = 2:end_index-1
        SD = default_struct; SD.subs{i_objdim} = i;
        SL = default_struct; SL.subs{i_objdim} = i-1;
        SH = default_struct; SH.subs{i_objdim} = i+1;
        der = (squeeze(subsref(obj.values,SH)) * (pts(i) - pts(i-1))^2 + ...
               squeeze(subsref(obj.values,SD)) * ((pts(i+1) - pts(i))^2 - (pts(i) - pts(i-1))^2) - ...
               squeeze(subsref(obj.values,SL)) * (pts(i+1) - pts(i))^2) / ...
              ((pts(i+1) - pts(i-1)) * (pts(i+1) - pts(i)) * (pts(i) - pts(i-1)));
        dCdD = subsasgn(dCdD, SD, der);
    end
    
    % Création de l'instance de sortie
    out = MultiDimVar(dCdD, obj.dim_names, obj.dim_points);

end