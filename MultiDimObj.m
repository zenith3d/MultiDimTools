classdef MultiDimObj
    
    properties
        data       = {};	% Table contenant les données
        shape      = [];	% Taille du tableau
        n_dims 	   = 0;		% Nombre de dimensions du tableau
        dim_names  = {};	% Noms des interpolants associés à chaque dimension
        dim_points = {};	% Valeurs des interpolants associés à chaque dimension
    end
    
    
    methods
        
        function obj = MultiDimObj(data, dim_names, dim_points)
            obj.data 	   = data;
            obj.shape 	   = size(values);
            obj.n_dims 	   = ndims(values);
            obj.dim_names  = dim_names;
            obj.dim_points = dim_points;
        end
        
        
        % Fonctions élémentaires
        % ----------------------------------------------------------------
        
        
        
        % Fonctions avancées
        % ---------------------------------------------------------------------
        
        
        
        % Overriding des opérateurs élémentaires
        % ----------------------------------------------------------------
        
        
    end
    
end