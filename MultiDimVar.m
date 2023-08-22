classdef MultiDimVar
    
    properties
        values     = [];	% Table contenant les valeurs de la variable
        shape      = [];	% Taille du tableau
        n_dims 	   = 0;		% Nombre de dimensions du tableau
        dim_names  = {};	% Noms des interpolants associ�s � chaque dimension
        dim_points = {};	% Valeurs des interpolants associ�s � chaque dimension
    end
    
    
    methods
        
        function obj = MultiDimVar(values, dim_names, dim_points)
            
            % Reshaping des tableaux contenant les interpolants
            dim_points = cellfun(@(x) reshape(x,1,[]), dim_points, ...
                'UniformOutput', false);
                        
            % Initialisation des propri�t�s
            % Reshaping du tableau 'values' dans le cas o� il est de
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

        end
        
        
        % Fonctions �l�mentaires
        % ----------------------------------------------------------------
           
        function out = squeeze(obj)
            % SQUEEZE Supprime les dimensions singuli�res (singleton)
            
            % D�termination des dimensions singuli�res
            singleton_dims = find(obj.shape == 1);
            
            % Si aucune dimension n'est singuli�re, on retourne l'instance
            % intouch�e
            if isempty(singleton_dims)
                out = obj;
                return;
            end
            
            % Suppression des dimensions singuli�res
            squeezed_dims       = setdiff(1:obj.n_dims, singleton_dims);
            squeezed_dim_names  = obj.dim_names(squeezed_dims);
            squeezed_dim_points = obj.dim_points(squeezed_dims);
            
            % Reshape du tableau
            squeezed_shape  = obj.shape(squeezed_dims);
            squeezed_values = reshape(obj.values, squeezed_shape);
            
            % Cr�ation de l'instance retourn�e
            out = MultiDimVar(squeezed_values, squeezed_dim_names, ...
                      squeezed_dim_points);
        end
        
        function out = augment(obj, new_dim_names, new_dim_points)
            % AUGMENT Augmente les dimensions de la variable
            
            % V�rification de la coh�rence des arguments
            if length(new_dim_names) ~= length(new_dim_points)
                error(['Augmentation impossible. Les arguments ne ' ...
                       'sont pas coh�rents entre eux : tailles ' ...
                       'diff�rentes']);
            end
            
            % D�termination des nouvelles dimensions
            [aug_dim_names, i_objdims, i_dims] = ...
                union(obj.dim_names, new_dim_names, 'stable');
            aug_dim_points = {obj.dim_points{i_objdims}, ...
                new_dim_points{i_dims}};
            aug_shape 	   = cellfun(@length, aug_dim_points);              % A REMPLACER PAR [obj.shape(i_objdims) ????]
            
            % Augmentation du tableau gr�ce � 'repmat'
            r = aug_shape; r(1:obj.n_dims) = 1;
            aug_values = repmat(obj.values, r);
            
            % Cr�ation de l'instance retourn�e
            out = MultiDimVar(aug_values, aug_dim_names, aug_dim_points);
        end
        
        function out = permute(obj, permute_dim_names)
            % PERMUTE Permute les dimensions de la variable de mani�re � avoir
            % ses dimensions dans l'ordre d�sir�
            
            % R�cup�ration des noms de dimensions qui concernent la
            % variable
            reord_dim_names = ...
                intersect(permute_dim_names, obj.dim_names, 'stable');
            if length(reord_dim_names) ~= length(obj.dim_names)
                error(['Permutation des dimensions impossible. La variable ' ...
                       'concern�e dispose d''une ou de plusieurs ' ...
                       'dimension(s) non contenue(s) dans le tableau ' ...
                       'ordonn�.']);
            end
            
            % D�termination des indices de permutation puis permutation des
            % noms et des valeurs associ�es
            [~, i_dims] = ismember(reord_dim_names, obj.dim_names);
            reord_dim_names  = obj.dim_names(i_dims);
            reord_dim_points = obj.dim_points(i_dims);
            
            % Permutation du tableau
            reord_values = permute(obj.values, i_dims);
            
            % Cr�ation de l'instance retourn�e
            out = MultiDimVar(reord_values, reord_dim_names, reord_dim_points);
        end
        
        function out = sortdims(obj)
            % SORTDIMS Trie les dimensions selon leur nom de mani�re croissante

            % Tri des dimensions selon leur nom
            [sorted_dim_names, i_dims] = sort(obj.dim_names);
            sorted_dim_points = obj.dim_points(i_dims);

            % Permutation du tabeau gr�ce � 'permute'
            sorted_values = permute(obj.values, i_dims);

            % Cr�ation de l'isntance retourn�e
            out = MultiDimVar(sorted_values, sorted_dim_names, sorted_dim_points);
        end
        
        function out = augmentas(obj, v)
            % AUGMENTAS Augmente les dimensions du tableau de mani�re � avoir
            % ses dimensions identiques � la variable 'v'
            out = obj.augment(v.dim_names, v.dim_points);
        end
        
        function out = permuteas(obj, v)
            % PERMUTEAS Permute les dimensions de la variable de mani�re � avoir
            % ses dimensions dans l'ordre de celles de la variable d�sir�e
            
            if obj.n_dims ~= v.n_dims || obj.n_dims == 1
                out = obj;
            else
                out = obj.permute(v.dim_names);
            end
        end
        
        
        % Fonctions avanc�es
        % ---------------------------------------------------------------------
        
        function out = interpn(obj, interp_dim_names, ...
                interp_dim_points, method, opts)
            % INTERPN Interpolation de la variable selon la ou les dimensions 
            % d�sir�es. Les tableaux contenus dans le cellarray 
            % 'interp_dim_points' doivent �tre sous forme de vecteurs.
            
            % V�rification de la coh�rence des arguments
            n_names  = length(interp_dim_names);
            n_points = length(interp_dim_points);
            if n_names ~= n_points
                error(['Interpolation impossible. Les arguments ne ' ...
                       'sont pas coh�rents entre eux : tailles ' ...
                       'diff�rentes']);
            end
            
            % D�termination des dimensions concern�es par l'interpolation :
            % intersection des dimensions
            [~, i_objdims, i_intdims] = ...
                intersect(obj.dim_names, interp_dim_names, 'stable');
            
            % Si aucune dimension n'est concern�e par l'interpolation, on
            % retourne la variable intouch�e
            if isempty(i_objdims)
                out = obj;
                return;
            end
            
            % R�cup�ration des points associ�s � chaque dimension pour
            % l'interpolation
            int_dim_names  = obj.dim_names;
            int_dim_points = obj.dim_points;
            for dim = 1:length(i_objdims)
                int_dim_points{i_objdims(dim)} = ...
                    interp_dim_points{i_intdims(dim)};
            end
            
            % Interpolation de la variable selon la m�thode choisie
            v = cell(1,length(int_dim_points));
            [v{:}] = ndgrid(int_dim_points{:});
            int_values = interpn(obj.dim_points{:}, obj.values, v{:}, method);
            
            % Cr�ation de l'instance de sortie
            out = MultiDimVar(int_values, int_dim_names, int_dim_points);
            
            % Suppression des dimensions singuli�res si d�sir�
            if nargin > 3 && strcmp(opts, 'squeeze')
                out = out.squeeze();
            end
        end
        
        function out = interpt(obj, interp_dim_names, ...
                interp_dim_points, method, opts)
            % INTERPT Interpolation de la variable selon la ou les dimensions
            % d�sir�es. Les tableaux contenus dans le cellarray
            % 'interp_dim_points' doivent �tre sous la forme de vecteurs de
            % dimensions �gales (en fonction du temps par ex.)
            
            % V�rification de la coh�rence des arguments
            n_names  = length(interp_dim_names);
            n_points = length(interp_dim_points);
            if n_names ~= n_points
                error(['Interpolation impossible. Les arguments ne ' ...
                       'sont pas coh�rents entre eux : tailles ' ...
                       'diff�rentes']);
            end
            
            % D�termination des dimensions concern�es par l'interpolation :
            % intersection des dimensions
            [~, i_objdims, i_intdims] = ...
                intersect(obj.dim_names, interp_dim_names, 'stable');
            
            % Si aucune dimension n'est conc�rn�e par l'interpolation, on
            % retourne la variable intouch�e
            if isemtpy(i_objdims)
                out = obj;
                return;
            end
            
            % R�cup�ration des valeurs associ�es � chaque dimension pour
            % l'interpolation
            int_dim_names  = obj.dim_names;
            int_dim_points = obj.dim_points;
            for dim = 1:length(i_objdims)
                int_dim_points{i_objdims(dim)} = ...
                    interp_dim_points{i_intdims(dim)};
            end
            
            % Interpolation de la variable selon la m�thode choisie
            v = cell(1,length(int_dim_points));
            [v{:}] = ndgrid(int_dim_points{:});
            int_values = interpn(obj.dim_points{:}, obj.values, v{:}, method);
            
            % Cr�ation de l'instance de sortie
            out = MultiDimVar(int_values, int_dim_names, int_dim_points);
            
            % Suppression des dimensions singuli�res si d�sir�
            if nargin > 3 && strcmp(opts, 'squeeze')
                out = out.squeeze();
            end
        end
        
        function out = extrap(obj, extrap_dim_names, extrap_dim_points, ...
                opts)
            % EXTRAP Extrapolation de la variable selon la dimension d�sir�e.
            % Les tableaux contenus dans le cellarray 'extrap_dim_points' 
            % doivent �tre sous forme de vecteurs.
            out = [];
        end
        
        function out = diff(obj, diff_dim_name, fill_missing)
            % DIFF Calcule la d�riv�e de la variable en fonction de la
            % dimension d�sir�e
            
            % D�termination de la dimension concern�e
            [~, i_objdim, ~] = intersect(obj.dim_names, diff_dim_name, 'stable');
            
            % Si aucune dimension n'est conc�rn�e, la d�riv�e est nulle
            if isempty(i_objdim)
                out = obj;
                out.values(:) = 0;
                return;
            end
            
            % R�cup�ration des points correspondants � la dimension concern�e
            points    = obj.dim_points{i_objdim};
            end_index = obj.shape(i_objdim);
            
            % Initialisation de la table de sortie 
            dCdD = zeros(obj.shape);
            
            % D�finition de la structure 'subs' par d�faut
            default_struct.type = '()';
            default_struct.subs = repmat({':'}, 1, obj.n_dims);

            % D�riv�e en debut de table (ordre 1)
            SD = default_struct; SD.subs{i_objdim} = 1;
            SH = default_struct; SH.subs{i_objdim} = 2;
            der = (squeeze(subsref(obj.values,SH)) - ...
                   squeeze(subsref(obj.values,SD))) / ...
                  (points(2) - points(1));
            dCdD = subsasgn(dCdD, SD, der);

            % D�riv�e en fin de table (ordre 1)
            SD = default_struct; SD.subs{i_objdim} = end_index;
            SL = default_struct; SL.subs{i_objdim} = end_index-1;
            der = (squeeze(subsref(obj.values,SD)) - ...
                   squeeze(subsref(obj.values,SL))) / ...
                  (points(end_index) - points(end_index-1));
            dCdD = subsasgn(dCdD, SD, der);

            % D�riv�e sur le reste de la table (ordre 2 - diff�rence centrale pond�r�e)
            for i = 2:end_index-1
                SD = default_struct; SD.subs{i_objdim} = i;
                SL = default_struct; SL.subs{i_objdim} = i-1;
                SH = default_struct; SH.subs{i_objdim} = i+1;
                der = (squeeze(subsref(obj.values,SH)) * (points(i) - points(i-1)) + ...
                       squeeze(subsref(obj.values,SD)) * (points(i+1) + points(i-1) - 2*points(i)) + ...
                       squeeze(subsref(obj.values,SL)) * (points(i) - points(i+1))) / ...
                      (2 * (points(i+1) - points(i)) * (points(i) - points(i-1)));
                dCdD = subsasgn(dCdD, SD, der);
            end
            
            % Cr�ation de l'instance de sortie
            out = MultiDimVar(dCdD, obj.dim_names, obj.dim_points);

        end
        
        function out = fillmissing(obj, dim_name, method)
            % FILLMISSING Remplie les entr�es manquantes dans le tableau
            
            % D�termination de la dimension concern�e
            [~, i_objdim, ~] = intersect(obj.dim_names, dim_name, 'stable');
            
            % V�rification
            if isempty(i_objdim)
                out = obj;
                return;
            end
            
            % Remplissage des valeurs manquantes selon la m�thode choisie
            out = obj;
            out.values = fillmissing(obj.values, method, i_objdim, ...
                                     'EndValues', 'extrap', ...
                                     'SamplePoints', obj.dim_points{i_objdim});
        end
        
        function [out, I] = minalongdims(obj, dim_names)
            % MINALONGDIMS D�termine les �l�ments minimums sur les dimensions 
            % sp�cifi�es par 'dim_names'
                        
            % D�termination des dimensions concern�es
            [~, i_dims, ~] = intersect(obj.dim_names, dim_names, 'stable');
            
            % V�rification
            if isempty(i_dims)
                out = [];
                return;
            end
                      
            % R�cup�ration des dimensions non concern�es
            i_not_dims     = setxor(1:obj.n_dims, i_dims);
            not_dim_names  = obj.dim_names(i_not_dims);
            not_dim_points = obj.dim_points(i_not_dims);
            
            % D�termination des �l�ments minimums sur les dimensions
            % d�sir�es
            [~,I] = min(obj.values, [], i_dims, 'linear');
            indices = cell(obj.n_dims, 1);
            [indices{:}] = ind2sub(obj.shape, squeeze(I));
            
            % Allocation du nombre d'instances 'MultiDimVar' � retourner
            n   = length(i_dims);
            out = cell(n,1);
                        
            % R�cup�ration des valeurs de chaque dimension qui minimise
            % 'obj.values'
            for i = 1:n
                i_dim  = i_dims(i);
                min_index_values = obj.dim_points{i_dim}(indices{i_dim});
                
                % Cr�ation de nouvelles instances MultiDimVar
                out{i} = MultiDimVar(min_index_values, not_dim_names, ...
                    not_dim_points);
            end
            
        end
        
        
        
        % Overriding des op�rateurs �l�mentaires
        % ----------------------------------------------------------------
        
        function r = lt(a,b)
            % LT (override) Less than (a < b)
            if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
                r = (a < b.values);
                return;
            elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
                r = (a.values < b);
                return;
            end
        end
        
        function r = uplus(v)
            % UPLUS (override) Unitary plus (+a)
            r = v;
        end
        
        function r = uminus(v)
            % UMINUS (override) Unitary minus (-a)
            r = v;
            r.values = -r.values;
        end
        
        function r = plus(a, b)
            % PLUS (override) Binary addition (a+b)
            
            %
            if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
                r = b;
                r.values = a + b.values;
                return;
            elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
                r = a;
                r.values = a.values + b;
                return;
            end
                        
            % Augmentation des variables si n�cessaire
            a_aug = a.augmentas(b);
            b_aug = b.augmentas(a);
            
            % Permutation des variables de l'une des deux variables
            b_aug = b_aug.permuteas(a_aug);
            
            % Addition des valeurs
            r = a_aug;
            r.values = a_aug.values + b_aug.values;
        end
        
        function r = minus(a, b)
            % MINUS (override) Binary substraction (a-b)
            r = a + uminus(b);
        end
        
        function r = times(a, b)
            % TIMES (override) Element-wise multiplication (a.*b)
            
            if ~isa(a,'MultiDimVar') || ~isa(b,'MultiDimVar')
                r = a * b;
                return;
            end
                        
            % Augmentation des variables si n�cessaire
            a_aug = a.augmentas(b);
            b_aug = b.augmentas(a);
            
            % Permutation des variables de l'une des deux variables
            b_aug = b_aug.permuteas(a_aug);
            
            % Addition des valeurs
            r = a_aug;
            r.values = a_aug.values .* b_aug.values;
        end
        
        function r = mtimes(a, b)
            % MTIMES (override) Matrix multiplication (a*b)
            
            % Multiplication par un scalaire
            if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
                r = b;
                r.values = a * b.values;
            elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
                r = a;
                r.values = a.values * b;
            else
                error('Type de multiplication non pris en charge');
            end
        end
        
        function r = rdivide(a, b)
            % RDIVIDE (override) Right element-wise division (a./b)
            
            if ~isa(a,'MultiDimVar') || ~isa(b,'MultiDimVar')
                r = a / b;
                return;
            end
                        
            % Augmentation des variables si n�cessaire
            a_aug = a.augmentas(b);
            b_aug = b.augmentas(a);
            
            % Permutation des variables de l'une des deux variables
            b_aug = b_aug.permuteas(a_aug);
            
            % Addition des valeurs
            r = a_aug;
            r.values = a_aug.values ./ b_aug.values;
        end
        
        function r = mrdivide(a, b)
            % MRDIVIDE (override) Matrix right division (a/b)
            
            % Division par un scalaire
            if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
                r = b;
                r.values = a ./ b.values;
            elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
                r = a;
                r.values = a.values / b;
                
            else
                error('Type de division non pris en charge');
            end
        end
        
        function r = power(a, b)
            % POWER (override) Element-wise power (a.^b)
            
            if ~isa(a,'MultiDimVar') || ~isa(b,'MultiDimVar')
                r = a ^ b;
                return;
            end
            
            % V�rification de la coh�rence entre les variables
%             if isempty(setxor(a.dim_names, b.dim_names))
%                 r = a;
%                 r.values = a.values .^ b.values;
%                 return;
%             end
            
            % Augmentation des variables si n�cessaire
            a_aug = a.augmentas(b);
            b_aug = b.augmentas(a);
            
            % Permutation des variables de l'une des deux variables
            b_aug = b_aug.permuteas(a_aug);
            
            % Addition des valeurs
            r = a_aug;
            r.values = a_aug.values .^ b_aug.values;
        end
        
        function r = mpower(a, b)
            % MRDIVIDE (override) Matrix power (a^b)
            
            % Multiplication par un scalaire
            if isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
                r = a;
                r.values = a.values .^ b;
            else
                error('Type de puissance non pris en charge');
            end
        end
        
        function r = cos(v)
            % COS (override) Cosinus (cos(v))
            r = v;
            r.values = cos(v.values);
        end
        
        function r = acos(v)
            % ACOS (override) Arccosinus (acos(v))
            r = v;
            r.values = acos(v.values);
        end
        
        function r = sin(v)
            % SIN (override) Sinus (sin(v))
            r = v;
            r.values = sin(v.values);
        end
        
        function r = asin(v)
            % ASIN (override) Arcsinus (asin(v))
            r = v;
            r.values = asin(v.values);
        end
        
        function r = tan(v)
            % TAN (override) Tangent (tan(v))
            r = v;
            r.values = tan(v.values);
        end
        
        function r = atan(v)
            % ATAN (override) Arctangent (atan(v))
            r = v;
            r.values = atan(v.values);
        end
        
        function r = atan2(a, b)
            % ATAN2 (override) Arctangent2 (atan2(a,b))
            r = a;
            r.values = atan2(a.values, b.values);
        end
        
        function r = sqrt(v)
            % SQRT (override) Square-root (sqrt(v))
            r = v;
            r.values = sqrt(v.values);
        end
        
        function r = max(a, b)
            % MAX (override) Maximum (max(a, b))
            if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
                r = b;
                r.values = max(a, b.values);
                return;
            elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
                r = a;
                r.values = max(a.values, b);
                return;
            end
        end
        
        function r = min(a, b)
            % MIN (override) Minimum (min(a, b))
            if ~isa(a,'MultiDimVar') && isa(b,'MultiDimVar')
                r = b;
                r.values = min(a, b.values);
                return;
            elseif isa(a,'MultiDimVar') && ~isa(b,'MultiDimVar')
                r = a;
                r.values = min(a.values, b);
                return;
            end
        end
        
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
    
end


% _____________________________________________________________________________
%
%                               SUBFUNCTIONS
% _____________________________________________________________________________

function err = inputs_consistency(values, dim_names, dim_points)

% V�rification de la coh�rence des entr�es
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
    err = 0;
end

end