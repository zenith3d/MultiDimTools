function out = solvealongdim(obj, dim_name, cost_var)
    % SOLVEALONGDIM Détermine la valeur de l'interpolant qui annule
    % la variable tout en minimisant la fonction de coût annexe
    % 'cost_var'.
    
    % Détermination de la dimension concernée par la résolution
    [~, i_objdim, ~] = ...
        intersect(obj.dim_names, dim_name, 'stable');
    
    % Si aucune dimension n'est concernée par l'interpolation, on
    % ne retourne rien
    if isempty(i_objdim)
        out = [];
        return;
    end
    
    %
    error('Fonction non implémentée');
end