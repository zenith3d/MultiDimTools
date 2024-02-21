function out = extract(obj, dim_names, dim_points, method)
    % EXTRACT Réalise l'extraction (interpolation / extrapolation)
    % de la variable selon la ou les dimensions spécifiées dans le 
    % cellarray 'dim_names'. Les tableaux contenus dans le 
    % cellarray 'dim_points' doivent être sous forme de vecteurs.
    
    % Vérification de la cohérence des arguments
    n_names  = length(dim_names);
    n_points = length(dim_points);
    if n_names ~= n_points
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
    
    % Récupération des points associés à chaque dimension pour
    % l'interpolation
    int_dim_names  = obj.dim_names;
    int_dim_points = obj.dim_points;
    for dim = 1:length(i_objdims)
        int_dim_points{i_objdims(dim)} = ...
            dim_points{i_intdims(dim)};
    end
    
    % Interpolation de la variable selon la méthode choisie
    xq = cell(1,length(int_dim_points));
    [xq{:}] = ndgrid(int_dim_points{:});
    int_values = interpn(obj.dim_points{:}, obj.values, xq{:}, ...
        method);                                                    % A AMELIORER : Utilisation directe de griddedInterpolant (voir interpn.m) (permet également d'extrapoler...)
    
    % Création de l'instance de sortie et suppression des
    % dimensions singulières
    out = MultiDimVar(int_values, int_dim_names, int_dim_points);
    out = out.squeeze();
end