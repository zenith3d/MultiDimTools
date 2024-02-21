function [out, I] = minalongdims(obj, dim_names, nanflag)
    % MINALONGDIMS Détermine les éléments minimums sur les 
    % dimensions spécifiées par 'dim_names'.
    
    if nargin < 3
        nanflag = 'omitnan';
    end
    
    % Détermination des dimensions concernées
    [~, i_objdims, ~] = ...
        intersect(obj.dim_names, dim_names, 'stable');
    
    % Vérification
    if isempty(i_objdims)
        out = {};
        I   = [];
        return;
    end
    
    % Récupération des dimensions non concernées
    i_not_dims     = setxor(1:obj.n_dims, i_objdims);
    not_dim_names  = obj.dim_names(i_not_dims);
    not_dim_points = obj.dim_points(i_not_dims);
    
    % Détermination des éléments minimums sur les dimensions
    % désirées
    [~,I] = min(obj.values, [], i_objdims, nanflag, 'linear');
    indices = cell(obj.n_dims, 1);
    [indices{:}] = ind2sub(obj.shape, squeeze(I));
    
    % Allocation du nombre d'instances 'MultiDimVar' à retourner
    n   = length(i_objdims);
    out = cell(n,1);
    
    % Récupération des valeurs de chaque dimension qui minimise
    % 'obj.values'
    for i = 1:n
        i_dim  = i_objdims(i);
        min_index_values = obj.dim_points{i_dim}(indices{i_dim});
        
        % Création de nouvelles instances MultiDimVar
        out{i} = MultiDimVar(min_index_values, not_dim_names, ...
            not_dim_points);
    end
    
end