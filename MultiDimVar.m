classdef MultiDimVar
    
    properties
        values     = [];	% Table contenant les valeurs de la variable
        shape      = [];	% Taille du tableau
        n_dims 	   = 0;		% Nombre de dimensions du tableau
        dim_names  = {};	% Noms des interpolants associés à chaque dimension
        dim_points = {};	% Valeurs des interpolants associés à chaque dimension
    end
    
    
    methods
        
        function obj = MultiDimVar(values, dim_names, dim_points, nochk)
            
            % Reshaping des tableaux contenant les interpolants
            dim_points = cellfun(@(x) reshape(x,[],1), dim_points, ...
                'UniformOutput', false);
                        
            % Initialisation des propriétés
            % Reshaping du tableau 'values' dans le cas où il est de
            % dimension 1
            if length(dim_names) == 1
                obj.values     = reshape(values,[],1);
                obj.shape 	   = length(values);
                obj.n_dims 	   = 1;
                obj.dim_names  = dim_names;
                obj.dim_points = dim_points;
            else
                obj.values 	   = values;
                obj.shape 	   = size(values);
                obj.n_dims 	   = ndims(values);
                obj.dim_names  = dim_names;
                obj.dim_points = dim_points;
            end

            % Vérification de la cohérence des entrées dans le cas d'une
            % création manuelle
            if nargin < 4 || ~nochk
                if (obj.n_dims ~= numel(obj.dim_names)) || ...
                   (obj.n_dims ~= numel(obj.dim_points)) || ...
                   any(obj.shape ~= cellfun(@numel, obj.dim_points))
                    error(['Impossible de créer une instance de ', ...
                    'MultiDimVar : dimensions incohérentes.']);
                end
            end

        end
        
        
        % Fonctions élémentaires
        % ----------------------------------------------------------------
        
        function out = squeeze(obj)
            % SQUEEZE Supprime les dimensions singulières (singletons)
            
            % Détermination des dimensions singulières
            singleton_dims = find(obj.shape == 1);
            
            % Si aucune dimension n'est singulière, on retourne une copie 
            % de l'instance intouchée
            if isempty(singleton_dims)
                out = obj;
                return;
            end
            
            % Suppression des dimensions singulières
            squeezed_dims       = setdiff(1:obj.n_dims, singleton_dims);
            squeezed_dim_names  = obj.dim_names(squeezed_dims);
            squeezed_dim_points = obj.dim_points(squeezed_dims);
            
            % Reshape du tableau
            squeezed_shape  = obj.shape(squeezed_dims);
            squeezed_values = reshape(obj.values, squeezed_shape);
            
            % Création de l'instance retournée
            out = MultiDimVar(squeezed_values, squeezed_dim_names, ...
                squeezed_dim_points, true);

        end
        
        function out = augment(obj, new_dim_names, new_dim_points)
            % AUGMENT Ajoute des nouvelles dimensions à la variable
            
            % Vérification de la cohérence des arguments
            if numel(new_dim_names) ~= numel(new_dim_points)
                error(['Augmentation impossible. Les arguments ne ' ...
                       'sont pas cohérents entre eux : tailles ' ...
                       'différentes.']);
            end
            
            % Détermination des nouvelles dimensions
            [aug_dim_names, i_objdims, i_newdims] = ...
                union(obj.dim_names, new_dim_names, 'stable');
            aug_dim_points = {obj.dim_points{i_objdims}, ...
                new_dim_points{i_newdims}};
            aug_shape 	   = cellfun(@length, aug_dim_points);              % A REMPLACER PAR [obj.shape(i_objdims) ????]
            
            % Augmentation du tableau grâce à 'repmat'
            r = aug_shape; r(1:obj.n_dims) = 1;
            aug_values = repmat(obj.values, r);
            
            % Création de l'instance retournée
            out = MultiDimVar(aug_values, aug_dim_names, ...
                aug_dim_points, true);
            
        end
        
        function out = permute(obj, perm_dim_names)
            % PERMUTE Permute les dimensions de la variable de manière à 
            % avoir ses dimensions dans l'ordre désiré
            
            % Récupération des noms de dimensions qui concernent la
            % variable
            [reord_dim_names, ~, i_objdims] = ...
                intersect(perm_dim_names, obj.dim_names, 'stable');
            if length(i_objdims) ~= obj.n_dims
                error(['Permutation des dimensions impossible. La ' ...
                       'variable concernée dispose d''une ou de ' ...
                       'plusieurs dimension(s) non contenue(s) dans ' ...
                       'le tableau ordonné fourni.']);
            end
            
            % Permutation inutile si une seule dimension
            if obj.n_dims < 2
                out = obj;
                return;
            end

            % Permutation des interpolants et du tableau
            reord_dim_points = obj.dim_points(i_objdims);
            reord_values = permute(obj.values, i_objdims);
            
            % Création de l'instance retournée
            out = MultiDimVar(reord_values, reord_dim_names, ...
                reord_dim_points, true);
            
        end
        
        function out = sortdims(obj)
            % SORTDIMS Trie les dimensions selon leur nom de manière 
            % croissante

            % Tri inutile si une seule dimension
            if obj.n_dims < 2
                out = obj;
                return;
            end

            % Tri des dimensions selon leur nom
            [sorted_dim_names, i_dims] = sort(obj.dim_names);
            sorted_dim_points = obj.dim_points(i_dims);

            % Permutation du tabeau grâce à 'permute()'
            sorted_values = permute(obj.values, i_dims);

            % Création de l'isntance retournée
            out = MultiDimVar(sorted_values, sorted_dim_names, ...
                sorted_dim_points, true);

        end
        
        function out = augmentas(obj, v)
            % AUGMENTAS Augmente les dimensions du tableau de manière à 
            % avoir ses dimensions identiques à la variable 'v' fournie
            
            out = obj.augment(v.dim_names, v.dim_points);
        end
        
        function out = permuteas(obj, v)
            % PERMUTEAS Permute les dimensions de la variable de manière à 
            % avoir ses dimensions dans un ordre identique à celles de la 
            % variable 'v' fournie
            
            out = obj.permute(v.dim_names);
        end
        
        function out = any(obj, dim_names)
            % ANYZERO Détermine si au moins un élément de la table est
            % non-nul
            
            if nargin < 2
                out = any(obj.values, 'all');
                return;
            end
            error('Fonction non implémentée');
        end
        
        function out = anynan(obj, dim_names)
            % ANYNAN Détermine si au moins un élément de la table est NaN
            
            if nargin < 2
                out = any(isnan(obj.values), 'all');
                return;
            end
            error('Fonction non implémentée');
        end
        
        
        % Fonctions avancées   
        % ----------------------------------------------------------------
              
        out = fillmissing(obj, dim_name, method)
        out = differentiate(obj, dim_name)
        out = integrate(obj, dim_name)
        out = extract(obj, dim_names, dim_points, method)
        out = interpn(obj, dim_names, dim_vars, method)
        out = solvealongdim(obj, dim_name, cost_var)
        [out, I] = minalongdims(obj, dim_names, nanflag)
        
        
        % Overriding des opérateurs élémentaires
        % ----------------------------------------------------------------
        
        function r = uplus(obj)
            % UPLUS (override) Unitary plus (+v)
            r = obj;
        end
        
        function r = uminus(obj)
            % UMINUS (override) Unitary minus (-v)
            r = obj;
            r.values = -r.values;
        end
        
        function r = ceil(obj)
            % CEIL (override) Round toward positive infinity
            r = obj;
            r.values = ceil(obj.values);
        end
        
        function r = floor(obj)
            % FLOOR (override) Round toward negative infinity
            r = obj;
            r.values = floor(obj.values);
        end
        
        function r = round(obj)
            % ROUND (override) Round to nearest integer
            r = obj;
            r.values = round(obj.values);
        end
        
        function r = cos(obj)
            % COS (override) Cosinus
            r = obj;
            r.values = cos(obj.values);
        end
        
        function r = acos(obj)
            % ACOS (override) Arccosinus
            r = obj;
            r.values = acos(obj.values);
        end
        
        function r = sin(obj)
            % SIN (override) Sinus
            r = obj;
            r.values = sin(obj.values);
        end
        
        function r = asin(obj)
            % ASIN (override) Arcsinus
            r = obj;
            r.values = asin(obj.values);
        end
        
        function r = tan(obj)
            % TAN (override) Tangent
            r = obj;
            r.values = tan(obj.values);
        end
        
        function r = atan(obj)
            % ATAN (override) Arctangent
            r = obj;
            r.values = atan(obj.values);
        end
        
        function r = exp(obj)
            % EXP (override) Exponential
            r = obj;
            r.values = exp(obj.values);
        end
        
        function r = log(obj)
            % LOG (override) Natural logarithm
            r = obj;
            r.values = log(obj.values);
        end
        
        function r = log10(obj)
            % LOG10 (override) Common logarithm (base 10)
            r = obj;
            r.values = log10(obj.values);
        end
        
        function r = sqrt(obj)
            % SQRT (override) Square-root
            r = obj;
            r.values = sqrt(obj.values);
        end
        
        function r = abs(obj)
            % ABS (override) Absolute value
            r = obj; 
            r.values = abs(obj.values);
        end
        
        function r = sign(obj)
            % SIGN (override) Sign function
            r = obj; 
            r.values = sign(obj.values);
        end
        
        function r = deg2rad(obj)
            % DEG2RAD (override) Convert angle from degrees to radians
            r = obj;
            r.values = deg2rad(obj.values);
        end
        
        function r = rad2deg(obj)
            % RAD2DEG (override) Convert angle from radians to degrees
            r = obj;
            r.values = rad2deg(obj.values);
        end
        
        function r = isnan(obj)
            % ISNAN (override) True for Not-a-Number
            r = obj;
            r.values = isnan(obj.values);
        end
        
        
        % Overriding des opérateurs élémentaires binaires (à deux entrées)
        % (https://fr.mathworks.com/help/matlab/ref/bsxfun.html)
        % ----------------------------------------------------------------
        
        function r = bsxfun(fun, a, b)
            % BSXFUN (override) Element-wise operation
            
            if isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
                
                % Augmentation des variables si nécessaire
                a = a.augmentas(b);
                b = b.augmentas(a);
                
                % Permutation des variables de l'une des deux variables
                b = b.permuteas(a);
                
                % Réalisation de l'opération
                r = a;
                r.values = fun(a.values, b.values);
                
            elseif ~isa(b,'MultiDimVar')
                r = a;
                r.values = fun(a.values, b);
            else
                r = b;
                r.values = fun(a, b.values);
            end
            
        end
        
        function r = plus(a, b)
            % PLUS (override) Binary addition (a + b)
            r = bsxfun(@plus, a, b);
        end
        
        function r = minus(a, b)
            % MINUS (override) Binary substraction (a - b)
            r = bsxfun(@minus, a, b);
        end
        
        function r = times(a, b)
            % TIMES (override) Element-wise multiplication (a .* b)
            r = bsxfun(@times, a, b);
        end
        
        function r = rdivide(a, b)
            % RDIVIDE (override) Right element-wise division (a ./ b)
            r = bsxfun(@rdivide, a, b);
        end
        
        function r = ldivide(a, b)
            % LDIVIDE (override) Left element-wise division (a .\ b)
            r = bsxfun(@ldivide, a, b);
        end
        
        function r = power(a, b)
            % POWER (override) Element-wise power (a .^ b)
            r = bsxfun(@power, a, b);
        end
                
        function r = max(a, b)
            % MAX (override) Maximum (max(a, b))
            r = bsxfun(@max, a, b);
        end
        
        function r = min(a, b)
            % MIN (override) Minimum (min(a, b))
            r = bsxfun(@min, a, b);
        end
        
        function r = rem(a, b)
            % REM (override) Remainder after division (rem(a, b))
            r = bsxfun(@rem, a, b);
        end
        
        function r = mod(a, b)
            % MOD (override) Modulus after division (mod(a, b))
            r = bsxfun(@mod, a, b);
        end
        
        function r = atan2(a, b)
            % ATAN2 (override) Four-quadrant inverse tangent (atan2(a,b))
            r = bsxfun(@atan2, a, b);
        end
        
        function r = eq(a, b)
            % EQ (override) Equal (a == b)
            r = bsxfun(@eq, a, b);
        end
        
        function r = ne(a, b)
            % NE (override) Not equal (a ~= b)
            r = bsxfun(@neq, a, b);
        end
        
        function r = lt(a, b)
            % LT (override) Less than (a < b)
            r = bsxfun(@lt, a, b);
        end
                
        function r = le(a, b)
            % LE (override) Less than or equal to (a <= b)
            r = bsxfun(@le, a, b);
        end
        
        function r = gt(a, b)
            % GT (override) Greater than (a > b)
            r = bsxfun(@gt, a, b);
        end
        
        function r = ge(a, b)
            % GE (override) Greater than or equal to (a >= b)
            r = bsxfun(@ge, a, b);
        end
        
        function r = and(a, b)
            % AND (override) Logical AND (A & B)
            r = bsxfun(@and, a, b);            
        end
        
        function r = or(a, b)
            % OR (override) Logical OR (A | B)
            r = bsxfun(@or, a, b);            
        end
        
        function r = xor(a, b)
            % OR (override) Logical OR (A | B)
            r = bsxfun(@xor, a, b);            
        end
        
        
            
%         function r = mtimes(a, b)
%             % MTIMES (override) Matrix multiplication (a*b)
%             
%             % Multiplication par un scalaire
%             if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
%                 r = b;
%                 r.values = a * b.values;
%             elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
%                 r = a;
%                 r.values = a.values * b;
%             else
%                 error('Type de multiplication non pris en charge');
%             end
%         end
        
%         function r = mrdivide(a, b)
%             % MRDIVIDE (override) Matrix right division (a/b)
%             
%             % Division par un scalaire
%             if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
%                 r = b;
%                 r.values = a ./ b.values;
%             elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
%                 r = a;
%                 r.values = a.values / b;
%                 
%             else
%                 error('Type de division non pris en charge');
%             end
%         end

%         function r = mpower(a, b)
%             % MRDIVIDE (override) Matrix power (a^b)
%             
%             % Multiplication par un scalaire
%             if isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
%                 r = a;
%                 r.values = a.values .^ b;
%             else
%                 error('Type de puissance non pris en charge');
%             end
%         end
                
%         function [M,I] = min(v)
%             % MIN (override) Minimum
%             
%             [M, min_index] = min(v.values(:));
%             indices = cell(v.n_dims,1);
%             [indices{:}] = ind2sub(v.shape, min_index);
%             I = cell2mat(indices);
%             
%         end
        
%         function out = subsref(v, s)
%             % SUBSREF (override) Subscripted reference
%             out = [];
%         end
%         
%         function out = subsasgn(a, s, b)
%             % SUBSASGN (override) Subscripted assignment
%             out = [];
%         end
        
    end
    
    methods (Static)
        
        % Constructeurs élémentaires
        % ----------------------------------------------------------------
        
        function out = zeros(dim_names, dim_values)
            % ZEROS (override) Zeros array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = zeros(shape);
            out = MultiDimVar(values, dim_names, dim_values, false);
        end
        
        function out = ones(dim_names, dim_values)
            % ONES (override) Ones array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = ones(shape);
            out = MultiDimVar(values, dim_names, dim_values, false);
        end
        
        function out = false(dim_names, dim_values)
            % FALSE (override) False array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = false(shape);
            out = MultiDimVar(values, dim_names, dim_values, false);
        end
        
        function out = true(dim_names, dim_values)
            % TRUE (override) True array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = true(shape);
            out = MultiDimVar(values, dim_names, dim_values, false);
        end
        
        function out = NaN(dim_names, dim_values)
            % NAN (override) Not-a-Number array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = NaN(shape);
            out = MultiDimVar(values, dim_names, dim_values, false);
        end
                
    end
    
    methods (Access = protected)
        
    end
end
