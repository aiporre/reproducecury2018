
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo / tutoriel pour le programme de matching %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pour toute question, mon email est : joan_glaunes@yahoo.fr

% Dans cet exemple on va apparier deux ensembles composés de quatre paires
% source/cible de types surfaces, courbes, nuages de points et landmarks

% 1/ Dans une première variable "s", on rentre tout ce qui concerne la
% partie déformations.

% s.x doit être un tableau de taille 3*n contenant les coordonées 3D de tous les
% points source. Ici la surface source est formée
% de deux triangles (4 sommets), la courbe source de quatre points,
% le nuage source comporte 2 points et il y a 2 landmarks source, donc n=12

setenv('LD_LIBRARY_PATH','')

clear s
s.x = [0,1,1,0,1.5,1.6,1.8,2, 3,3.2,3.8, 4;
       0,0,1,1, 0 , 0 , 0 ,0, 0, 0 , 0 , 0 ;
       0,0,0,0, 0 , 0 , 0 ,0, 0, 0 , 0 , 0];
   
% s.sigmaV est l'échelle de déformation, intervenant dans le calcul du noyau. Il
% faut fixer une valeur correspondant à l'ordre de grandeur des coordonnées des
% points. Par exemple, si ces coordonnées vont de 0 à 250, et que les
% courbes à
% apparier ne sont pas trop isolées dans une région de l'espace, sigmaV = 20 est
% un bon choix. Ici les coordonnées sont de l'ordre de l'unité, on va
% choisir sigmaV = 0.5
s.sigmaV = .5;

% Autres paramêtres importants

s.gammaR = 0;        % poids du terme de régularité dans la fonctionnelle
% Avec gammaR = 0, on minimise uniquement le terme d'attache aux données,
% mais la déformation reste quand même régulière du fait de l'espace dans
% lequel s'effectue la minimisation

s.numbminims = 1; % signifie que l'on effectue une seule minimisation. Pour améliorer
% le matching, on peut fixer par exemple numbminims=3 et l'algo effectue
% alors plusieurs minimisations en diminuant à chaque fois la taille des
% noyaux d'appariement sigmaW, sigmaI (cf plus bas)

s.usefgt = 0; % =1 signifie que l'on utilise la Fast Gauss Transform pour les convolutions 
% (efficace à partir de 50 points environ).


% 2/ Dans une deuxième variable "target", on met tous ce qui concerne les
% cibles et les méthodes d'appariement utilisées

clear target

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Première cible de type surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{1}.method = 'surfcurr';

% indices des triangles de la surface SOURCE dans target{1}.vx
% Ces indices font référence à des colonnes de s.x
target{1}.vx = [1,1;
                2,3;
                3,4];
% coordonnées des sommets de la surface CIBLE dans target{1}.y
target{1}.y = [.5,1,0
                0,1,.5
               .5,1,.5];
% indices des triangles de la surface CIBLE dans target{1}.vy            
target{1}.vy = [1;
                2;
                3];
% échelle du noyau d'appariement de surfaces. Plus cette valeur est petite,
% plus l'appariement est précis. Un bon choix est une valeur
% de l'ordre de la distance entre l'objet source et cible. Cependant si on
% effectue plusieurs minimisations successives (variable s.numbminims), on
% peut fixer une valeur petite car l'algorithme effectue les premières étapes
% à des échelles supérieures
target{1}.sigmaW = .5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deuxième cible de type courbe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{2}.method = 'curvecurr';

% indices des segments de la courbe SOURCE dans target{2}.vx
% Ces indices font référence à des colonnes de s.x
target{2}.vx = [5,6,7;
                6,7,8];
% coordonnées des points de la courbe CIBLE dans target{2}.y
target{2}.y = [1.5, 2 ,2.5;
                0 , 0 ,0.5;
               0.5,0.5, 1 ];
% indices des segments de la courbe CIBLE dans target{2}.vy            
target{2}.vy = [1,2;
                2,3];
% échelle du noyau d'appariement de courbes
target{2}.sigmaW = .5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Troisième cible de type nuage de points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{3}.method = 'measures';

% indices des points SOURCE dans target{3}.vx
% Ces indices font référence à des colonnes de s.x
target{3}.vx = [9,10];
% coordonnées des points CIBLE dans target{3}.y
target{3}.y = [3 , 3.2,3.4;
              0.5, 0.7,0.5;
               0.5,0.5,0.5];
% échelle du noyau d'appariement de nuages de points
target{3}.sigmaI = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quatrième cible de type landmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target{4}.method = 'landmarks';

% indices des points SOURCE dans target{4}.vx
% Ces indices font référence à des colonnes de s.x
target{4}.vx = [11,12];
% coordonnées des points CIBLE dans target{3}.y
target{4}.y = [ 4 , 4;
               0.5,0.5;
               0.5,0.8];

s.optim_verbosemode = 0;

% On lance le programme:
s = matchCpp(s,target);

% La structure retournée "s" contient entre autres les variables "X" 
% (trajectoires de points template) et "mom" (vecteurs moments)
% qui paramètrent la déformation optimale. Aussi "distIdPhi" donne le coût
% de déformation obtenu (racine carrée de l'énergie)
disp(['coût de déformation D = D(id,phi) = ',num2str(s.distIdPhi)])
disp(' ')

% affichage du résultat
s.show = {'phi','y','xrig','x'};
s.showtraj = 1;    % affiche les trajectoires des points template (en bleu)
s.showmomtraj = 0; % affiche les vecteurs moments (flèches vertes)
s.showpoints = 1;
s.showlegend = 0;
clf
affiche(s);
view(-20,20)
axis equal

% Le programme flow.m permet ensuite de calculer l'image d'un ensemble de
% points 3D par la déformation optimale

disp('L''image des points')
V = rand(3,4)
disp('par la transformation est')
phiV = flowCpp(s,V)


