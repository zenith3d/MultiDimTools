function out = interpn(obj, dim_names, dim_vars, method)
    % INTERPN Réalise l'interpolation de la variable selon les
    % dimensions spécifiées dans 'dim_names'. À chaque dimension
    % spécifiée correspond une instance MultiDimVar contenue dans
    % le cellarray 'dim_vars'.
    
    % Vérification de la cohérence des arguments
    n_names = length(dim_names);
    n_vars  = length(dim_vars);
    if n_names ~= n_vars
        error(['Interpolation impossible. Les arguments ne ' ...
               'sont pas cohérents entre eux : tailles ' ...
               'différentes.']);
    end
    
    % Détermination des dimensions concernées par l'interpolation :
    % intersection des dimensions
    [~, i_objdims, i_intdims] = ...
        intersect(obj.dim_names, dim_names, 'stable');
    
    % Si aucune dimension n'est concernée par l'interpolation, on
    % retourne la variable intouchée
    if isempty(i_objdims)
        out = obj;
        return;
    end
    
    % Augmentation de la variable pour être cohérente des variables
    % fournies en entrée
    aug_obj = obj;
    for i = 1:length(i_intdims)
        aug_obj = aug_obj.augmentas(dim_vars{i_intdims(i)});
    end
    
    % Ensuite, on réalise une augmentation puis une permutation des 
    % variables d'entrée
    aug_dim_vars = cell(1,length(i_intdims));
    for i = 1:length(i_intdims)
        aug_dim_vars{i_intdims(i)} = ...
            dim_vars{i_intdims(i)}.augmentas(aug_obj).permuteas(aug_obj);
    end
    
    % Création des grilles nécessaires pour l'interpolation
    xq = cell(1,aug_obj.n_dims);
    [xq{:}] = ndgrid(aug_obj.dim_points{:});
    for i = 1:length(i_intdims)
        xq{i_objdims(i)} = aug_dim_vars{i_intdims(i)}.values;
    end
    
    error('Fonction non implémentée');
    
    % Interpolation de la variable selon la méthode choisie
    int_values = interpn(aug_obj.dim_points{:}, aug_obj.values, ...
        xq{:}, method);                                          	% A AMELIORER : Utilisation directe de griddedInterpolant (voir interpn.m) (permet également d'extrapoler...)
    
    % Création de l'instance de sortie et suppression des
    % dimensions singulières
    out = MultiDimVar(int_values, int_dim_names, aug_obj.dim_points);
    out = out.squeeze();
end