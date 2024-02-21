classdef MultiDimVar
    
    properties
        values     = [];	% Table contenant les valeurs de la variable
        shape      = [];	% Taille du tableau
        n_dims 	   = 0;		% Nombre de dimensions du tableau
        dim_names  = {};	% Noms des interpolants associés à chaque dimension
        dim_points = {};	% Valeurs des interpolants associés à chaque dimension
    end
    
    
    methods
        
        function obj = MultiDimVar(values, dim_names, dim_points)
            
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

        end
        
        
        % Fonctions élémentaires
        %  - squeeze(...)
        %  - augment(...)
        %  - permute(...)
        %  - sortdims(...)
        %  - augmentas(...)
        %  - permuteas(...)
        %  - any(...)
        %  - anynan(...)
        % ----------------------------------------------------------------
           
        function out = squeeze(obj)
            % SQUEEZE Supprime les dimensions singulières (singletons)
            
            % Détermination des dimensions singulières
            singleton_dims = find(obj.shape == 1);
            
            % Si aucune dimension n'est singulière, on retourne l'instance
            % intouchée
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
                      squeezed_dim_points);
        end
        
        function out = augment(obj, new_dim_names, new_dim_points)
            % AUGMENT Ajoute des nouvelles dimensions à la variable
            
            % Vérification de la cohérence des arguments
            if length(new_dim_names) ~= length(new_dim_points)
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
            out = MultiDimVar(aug_values, aug_dim_names, aug_dim_points);
        end
        
        function out = permute(obj, permute_dim_names)
            % PERMUTE Permute les dimensions de la variable de manière à 
            % avoir ses dimensions dans l'ordre désiré
            
            % Récupération des noms de dimensions qui concernent la
            % variable
            reord_dim_names = ...
                intersect(permute_dim_names, obj.dim_names, 'stable');
            if length(reord_dim_names) ~= length(obj.dim_names)
                error(['Permutation des dimensions impossible. La ' ...
                       'variable concernée dispose d''une ou de ' ...
                       'plusieurs dimension(s) non contenue(s) dans ' ...
                       'le tableau ordonné fourni.']);
            end
            
            % Détermination des indices de permutation puis permutation des
            % noms et des valeurs associées
            [~, i_dims] = ismember(reord_dim_names, obj.dim_names);
            reord_dim_names  = obj.dim_names(i_dims);
            reord_dim_points = obj.dim_points(i_dims);
            
            % Permutation du tableau
            reord_values = permute(obj.values, i_dims);
            
            % Création de l'instance retournée
            out = MultiDimVar(reord_values, reord_dim_names, ...
                reord_dim_points);
        end
        
        function out = sortdims(obj)
            % SORTDIMS Trie les dimensions selon leur nom de manière 
            % croissante

            % Tri des dimensions selon leur nom
            [sorted_dim_names, i_dims] = sort(obj.dim_names);
            sorted_dim_points = obj.dim_points(i_dims);

            % Permutation du tabeau grâce à 'permute()'
            sorted_values = permute(obj.values, i_dims);

            % Création de l'isntance retournée
            out = MultiDimVar(sorted_values, sorted_dim_names, ...
                sorted_dim_points);
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
            
            if obj.n_dims ~= v.n_dims || obj.n_dims == 1
                out = obj;
            else
                out = obj.permute(v.dim_names);
            end
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
        %  - fillmissing(...)
        %  - diff(...)
        %  - integr(...)
        %  - extract(...)           /!\ A AMELIORER /!\
        %  - interpn(...)           /!\ A IMPLEMENTER /!\
        %  - solvealongdim(...)     /!\ A IMPLEMETER /!\
        %  - solvealongdims(...)    /!\ A IMPLEMETER /!\
        %  - minalongdims(...)      
        % ----------------------------------------------------------------
              
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
        
%         function out = solvealongdims(obj, dim_names, cost_var)
%            % SOLVEALONGDIMS Détermine la valeur des interpolants selon les 
%            % dimensions désirées 'dim_names' qui annulent la variable tout 
%            % en minimisant la fonction de coût annexe 'cost_var'.
%            
%             % Détermination des dimensions concernées par la résolution
%             [~, i_objdims, ~] = ...
%                 intersect(obj.dim_names, dim_names, 'stable');
%             
%             % Si aucune dimension n'est concernée par l'interpolation, on
%             % ne retourne rien
%             if isempty(i_objdims)
%                 out = {};
%                 return;
%             end
%             
%             % Récupération des dimensions non concernées
%             i_not_dims     = setxor(1:obj.n_dims, i_objdims);
%             not_dim_names  = obj.dim_names(i_not_dims);
%             not_dim_points = obj.dim_points(i_not_dims);
%             
% %             % Détermination des éléments minimums sur les dimensions
% %             % désirées
% %             [~,I] = min(obj.values, [], i_dims, nanflag, 'linear');
% %             indices = cell(obj.n_dims, 1);
% %             [indices{:}] = ind2sub(obj.shape, squeeze(I));
%             
%             % Allocation du nombre d'instances 'MultiDimVar' à retourner
%             n   = length(i_objdims);
%             out = cell(n,1);
%             
% %             % Récupération des valeurs de chaque dimension qui minimise
% %             % 'obj.values'
% %             for i = 1:n
% %                 i_dim  = i_dims(i);
% %                 min_index_values = obj.dim_points{i_dim}(indices{i_dim});
% %                 
% %                 % Création de nouvelles instances MultiDimVar
% %                 out{i} = MultiDimVar(min_index_values, not_dim_names, ...
% %                     not_dim_points);
% %             end
%         end
%         
%         function out = interpt(obj, interp_dim_names, ...
%                 interp_dim_points, method, opts)
%             % INTERPT Réalise l'interpolation de la variable selon la ou 
%             % les dimensions désirées, qui évoluent toutes en fonction d'un
%             % paramètre commun (en fonction du temps par ex.). Les tableaux
%             % contenus dans le cellarray 'interp_dim_points' doivent être 
%             % sous la forme de vecteurs de dimensions égales.
%             
%             % Vérification de la cohérence des arguments
%             n_names  = length(interp_dim_names);
%             n_points = length(interp_dim_points);
%             if n_names ~= n_points
%                 error(['Interpolation impossible. Les arguments ne ' ...
%                        'sont pas cohérents entre eux : tailles ' ...
%                        'différentes.']);
%             end
%             
%             % Détermination des dimensions concernées par l'interpolation :
%             % intersection des dimensions
%             [~, i_objdims, i_intdims] = ...
%                 intersect(obj.dim_names, interp_dim_names, 'stable');
%             
%             % Si aucune dimension n'est concérnée par l'interpolation, on
%             % retourne la variable intouchée
%             if isemtpy(i_objdims)
%                 out = obj;
%                 return;
%             end
% 
%             error('Fonction non implémentée'); 
%             % - Revoir les deux paragraphes suivants : construction de v{:} à changer
%             % - Ajouter de nouveaux args en entrée de la fonction : le nom du nouvel 
%             %   interpolant (par ex. 'TIME') ainsi que ses valeurs ;
% 
%             % Récupération des valeurs associées à chaque dimension pour
%             % l'interpolation
%             int_dim_names  = obj.dim_names;
%             int_dim_points = obj.dim_points;
%             for dim = 1:length(i_objdims)
%                 int_dim_points{i_objdims(dim)} = ...
%                     interp_dim_points{i_intdims(dim)};
%             end
%             
%             % Interpolation de la variable selon la méthode choisie
%             v = cell(1,length(int_dim_points));
%             [v{:}] = ndgrid(int_dim_points{:});
%             int_values = interpn(obj.dim_points{:}, obj.values, v{:}, ...
%                 method);
%             
%             % Création de l'instance de sortie
%             out = MultiDimVar(int_values, int_dim_names, int_dim_points);
%             
%             % Suppression des dimensions singulières si désiré
%             if nargin > 3 && strcmp(opts, 'squeeze')
%                 out = out.squeeze();
%             end
%         end
%         
%         function out = extrap(obj, extrap_dim_names, extrap_dim_points, ...
%                 opts)
%             % EXTRAP Réalise l'extrapolation de la variable selon la 
%             % dimension désirée. Les tableaux contenus dans le cellarray 
%             % 'extrap_dim_points' doivent être sous forme de vecteurs.
%             out = [];
%             error('Extrapolation impossible. Fonction non implémentée.')
%         end
        
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
                a_aug = a.augmentas(b);
                b_aug = b.augmentas(a);
                
                % Permutation des variables de l'une des deux variables
                b_aug = b_aug.permuteas(a_aug);
                
                % Réalisation de l'opération
                r = a_aug;
                r.values = fun(a_aug.values, b_aug.values);
                
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
            out    = MultiDimVar(values, dim_names, dim_values);
        end
        
        function out = ones(dim_names, dim_values)
            % ONES (override) Ones array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = ones(shape);
            out = MultiDimVar(values, dim_names, dim_values);
        end
        
        function out = false(dim_names, dim_values)
            % FALSE (override) False array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = false(shape);
            out    = MultiDimVar(values, dim_names, dim_values);
        end
        
        function out = true(dim_names, dim_values)
            % TRUE (override) True array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = true(shape);
            out    = MultiDimVar(values, dim_names, dim_values);
        end
        
        function out = NaN(dim_names, dim_values)
            % NAN (override) Not-a-Number array
            if length(dim_names) == 1
                shape = [length(dim_values{1}), 1];
            else
                shape = cellfun(@length, dim_values);
            end
            values = NaN(shape);
            out    = MultiDimVar(values, dim_names, dim_values);
        end
                
    end
    
    methods (Access = protected)
        
    end
end


% _____________________________________________________________________________
%
%                               SUBFUNCTIONS
% _____________________________________________________________________________

function err = inputs_consistency(values, dim_names, dim_points)

% Vérification de la cohérence des entrées
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
