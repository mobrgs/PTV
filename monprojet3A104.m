% mon projet Matlab de detection / suivi de particules
clear % efface le contenu du Workspace
close all % ferme toutes les figures
clc % efface le contenu de la fenetre de travail (Command Window)
%==========================================================================
% Parametres utilisateur
Var1 = 100; % Nombre d'image à traiter
VarN = 'NomFigure'; % Titre de la figure des résultats enregistrée sous forme .png
% Nothing to modify below this line
%==========================================================================
% filename : structure contenant les noms des images a traiter
% pathname : chemin d'acces au repertoire contenant les images a traiter
% C : vecteur contenant les coordonnées colonne des pixels des maximums locaux (image amont)
% L : vecteur contenant les coordonnées ligne des pixels des maximums locaux (image amont)
% C2 : vecteur contenant les coordonnées colonne des pixels des maximums locaux (image avale)
% L2 : vecteur contenant les coordonnée ligne des pixels des maximums locaux (image avales)
% C2_zoom : vecteur contenant les coordonnées colonne des pixels des maximums locaux dans une partie zoomée de l'image avale
% C2_zoom : vecteur contenant les coordonnées ligne des pixels des maximums locaux dans une partie zoomée de l'image avale
% dilate_I : image amont dont les maximums locaux sont étendus
% dilate_I2 : image avale dont les maximums locaux sont étendus
% dilate_zoom_I2 : partie zoomée de l'image avale dont les maximums locaux sont étendus
% ext_zoom_I2 : zoom_I2 complété de pixels noirs sur se périphérie de taille "len"
% I1 : matrice contenant les intensités des pixels de l'image amont
% I2 : matrice contenant les intensité des pixels de l'image avale
% I1_ext : I1 étendue d'une taille "len" de lignes et colonnes nulles
% I2_ext : I2 étendue d'une taille "len" de lignes et colonnes nulles
% mat_corr : matrice contenant les valeurs de corrélation entre les particules de l'image amont (lignes) et avales (colonnes)
% mat_sol : matrice dont la première colonne prend le numéro de particule de l'image amont et la deuxième colonne les numéro de particule de l'image avale correspondante
% N : taille du zoom effectiué de part et d'autre de chaque pixel pris comme maximum local
% N2 : taille du zoom dans la zone de déplacement d'une particule d'une image à l'autre
% part_coord_x : coordonnées x des centres des particules (lignes : numéro particule ; colonne : numéro image)
% part_coord_y : coordonnées y des centres des particules (lignes : numéro particule ; colonne : numéro image)
% part_init : particule regardée pour la corrélation sur l'image amont
% part_test : particule regardée pour la corrélation sur l'image avale
% pos : vecteur contenant les position linéaire des maximums locaux d'une image
% R : coefficient de corrélation
% seuil : seuil sous lequel les maximums locaux ne sont pas détectés (image entière)
% seuil2 : seuil sous lequel les maximums locaux ne sont pas détectés (partie zoomée de l'image avale)
% sigma : écart type de la fonction gaussienne
% Xsol : coordonnées x des positions exactes des paticules sur l'image amont
% Ysol : coordonnées y des positions exactes des paticules sur l'image amont
% Xsol2 : coordonnées x des positions exactes des paticules sur l'image avale
% Ysol2 : coordonnées y des positions exactes des paticules sur l'image avale
% zoom_I2 : partie zoomée de l'image avale pour la récherche des particule correspondantes à celles présente dans l'image amont
% zoom : zoom autour d'une particule à la distance N du pixel central
% compteur : pour contenant le nombre de fois où une boucle est exécutée

[filename, pathname] = uigetfile({' * .tiff; * .png'}, 'Select ...images','MultiSelect','on');
if isequal(filename,0)
disp('User selected Cancel')
return
end
cd(pathname); % rend le répertoire d'acces des images actif

load('parameters.mat') % charge les paramètres du problème

%% Lecture de la première image

 sigma = 1.;
 I = imread(filename{1}); % conversion de l'image en matrice remplie de son intensité en niveau de gris par pixel
 
%% Boucle d'initialisation pour avoir les coordonnées des particules à traquer depuis l'image 1
 
     I = imread(filename{1}); % conversion de l'image 1 en matrice remplie de son intensité en niveau de gris
     I=im2double(I); % mise en double précision des valeurs contenues dans la matrice
     len = 5; 
     I_ext = wextend(2,'zpd',I,len); % étendue de la matrice avec des ligne et colonne nulles sur une taille "len"
     
     SE = strel('disk', 3);
     dilate_I=imdilate(I_ext,SE); % étendue des intensités des maximums locaux sur la taille estimée d'une particule

     seuil = multithresh(dilate_I); % estimation du seuil pour repérer les maximums locaux
     pos=find((I_ext)==(dilate_I) & (I_ext)>=seuil); % détermination des positions linéaire de ces maximums locaux
     [L,C]=ind2sub(size(I_ext), pos); % convertion de ces positions linéaires en coordonnées ligne, colonne

     part_coord_x = zeros(size(L,1),Var1); % création de la matrice contenant les coordonnées x des particules pour toutes les images
     part_coord_y = zeros(size(L,1),Var1); % création de la matrice contenant les coordonnées y des particules pour toutes les images

     N=2; 

     Xsol=zeros(size(C,1),1); % création du vecteur solution des coordonnées x des particules pour l'image 1
     Ysol=zeros(size(L,1),1); % création du vecteur solution des coordonnées x des particules pour l'image 1

     % boucle permettant de calculer la position exacte des centres des particules
     for i = 1:size(L,1);
        zoom=I_ext(L(i)-N:L(i)+N,C(i)-N:C(i)+N);
        vect_sol=[0 0 0];
        [Imax]=max(zoom(:));
        [solution]=lsqnonlin(@(vect_sol)ecart(vect_sol,zoom,sigma),[N, N, Imax]); % fonction de minimisation en concidérant une répartition gaussienne de l'intensité autour du centre
         Xsol(i)=solution(1)+C(i)-N-1;
         Ysol(i)=solution(2)+L(i)-N-1;
     end
     
     % boucle remplissant la première colonne de la matrice enregistrant les coordonnées des particules
     % ce vont être les particules traquées
     for i=1:size(part_coord_x,1);
        part_coord_x(i,1)=Xsol(i);
        part_coord_y(i,1)=Ysol(i);
     end

%% Boucle sur toutes les images

for n=1:Var1-1;

    % boucle transferant les coordonnées des maximums locaux précédemment enregistrées pour l'image avale comme les coordonnées de l'images amont
    if (n~=1);
        L=zeros(size(part_coord_x,1));
        C=zeros(size(part_coord_x,1));
        for i=1:size(mat_sol,1);
            if (mat_sol(i,2)~=0);
                L(i)=L2(mat_sol(i,2));
                C(i)=C2(mat_sol(i,2));
            end
        end
    end  
    
     I1 = imread(filename{n}); % conversion de l'image amont en matrice remplie de son intensité en niveau de gris
     I2 = imread(filename{n+1}); % conversion de l'image avale en matrice remplie de son intensité en niveau de gris

     I1=im2double(I1); % mise en double précision des valeurs contenues dans la matrice
     I2=im2double(I2); % mise en double précision des valeurs contenues dans la matrice

     I1_ext = wextend(2,'zpd',I1,len); % étendue de la matrice avec des ligne et colonne nulles sur une taille "len"
     I2_ext = wextend(2,'zpd',I2,len); % étendue de la matrice avec des ligne et colonne nulles sur une taille "len"

     SE = strel('disk', 3);
     dilate_I2=imdilate(I2_ext,SE); % étendue des intensités des maximums locaux sur la taille estimée d'une particule de l'image avale
     seuil = multithresh(dilate_I2); % estimation du seuil pour repérer les maximums locaux
     pos=find((I2_ext)==(dilate_I2) & (I2_ext)>=seuil); % détermination des positions linéaire de ces maximums locaux
     [L2,C2]=ind2sub(size(I2_ext), pos); % convertion de ces positions linéaires en coordonnées ligne, colonne

     mat_corr=zeros(size(L,1), size(L2,1)); % création de la matrice de corrélation

     N2=5;

    % boucle calculant la corrélation entre une particule de l'image amont et des particules de l'image avale contenunes dans une zone prédéfinies
    for j=1:size(L,1);
        if (L(j)~=0);
         zoom_I2=I2_ext(L(j)-N2:L(j)+N2,C(j)-N2:C(j)+N2); % définition de la zone de recherche sur l'image avale
         SE2 = strel('disk', 3);
         % recherche des maximums locaux dans cette zone prédéfinies
         dilate_zoom_I2=imdilate(zoom_I2,SE2);  % étendue des intensités des maximums locaux sur la taille estimée d'une particule
         seuil2 = multithresh(dilate_zoom_I2);  % estimation du seuil pour repérer les maximums locaux
         pos2=find((zoom_I2)==(dilate_zoom_I2) & (zoom_I2)>=seuil2);  % détermination des positions linéaire de ces maximums locaux
         [L2_zoom,C2_zoom]=ind2sub(size(zoom_I2), pos2); % convertion de ces positions linéaires en coordonnées ligne, colonne

         part_init=I1_ext(L(j)-N:L(j)+N,C(j)-N:C(j)+N); % découpage autour de la particule de l'image amont
         ext_zoom_I2=wextend(2,'zpd',zoom_I2,len);

             % calcul de la coorélation de toutes les particules présentes dans la zone prédéfinie précédemment avec la particule de l'image amont
             for k=1:size(C2_zoom,1)
                 part_test=ext_zoom_I2(L2_zoom(k)+len-N:L2_zoom(k)+len+N,C2_zoom(k)+len-N:C2_zoom(k)+len+N); % découpage autour de la particule à tester dans la zone prédéfinie de l'image avale
                 R=corr2(part_init,part_test); % calcul du coefficient de coréélation entre les deux particules
                 num=find((C2==C2_zoom(k)+C(j)-N2-1)&(L2==L2_zoom(k)+L(j)-N2-1)); % recherche du numéro de la particule de l'image avale
                 mat_corr(j,num)=R; % remplissage de la matrice de corrélation (ligne : numéro particule image amont, colonne : numéro particule image avale)
             end
        end
    end

    % création de la matrice de correspondance entre les particules de l'image amont et de l'image avale
    mat_sol=zeros(min(size(mat_corr)),2);

    % boucle remplissant cette matrice de correspondance
    for k=1:size(mat_sol,1);
        if (any(mat_corr(k,:))==1);
            a=find(mat_corr(k,:)==max(mat_corr(k,:))); % recherche de la corrélation maximum entre les particules des deux images
            if ((mat_corr(k,a))==max(mat_corr(:,a))); % vérification de l'unicité de la solution
                mat_sol(k,1)=k; % remplissage de la première colonne de la matrice de correspondance
                mat_sol(k,2)=a; % remplissage de la deuxième colonne de la matrice de correspondance
            end
        end
    end

     Xsol2=zeros(size(C2,1),1); % création du vecteur enregistrant les coordonnées x des centres des particules
     Ysol2=zeros(size(L2,1),1); % création du vecteur enregistrant les coordonnées y des centres des particules

    % boucle permettant de calculer la position exacte des centres des particules
    for i = 1:size(L2,1)
        zoom2=I2_ext(L2(i)-N:L2(i)+N,C2(i)-N:C2(i)+N);
        vect_sol=[0 0 0];
        [Imax]=max(zoom2(:));
        [solution]=lsqnonlin(@(vect_sol)ecart(vect_sol,zoom2,sigma),[N, N, Imax]);
         Xsol2(i)=solution(1)+C2(i)-N-1;
         Ysol2(i)=solution(2)+L2(i)-N-1;
    end

    % remplissage des matrices contenant toutes les solutions 
    for i=1:size(mat_sol,1);
        if (mat_sol(i,2)~=0); % ne rentre dans la boucle que s'il existe une correspondance entre la particule de l'image amont et une de l'image avale
            part_coord_x(i,n+1)=Xsol2(mat_sol(i,2)); % coordonnées x
            part_coord_y(i,n+1)=Ysol2(mat_sol(i,2)); % coordonnées y
        end
    end 
 end

 % save("Projet_donnees"); % sauvegarde le workspace

%% Ecriture du graphique présentant les résultats de la méthode PTV

compteur = 0.; 

% boucle comptant le nombre de particule visible sur toutes les images
for i=1:size(part_coord_x,1);
     if (all(part_coord_x(i,:))==1); % ne rentre dans la boucle que si aucune case n'est vide
         compteur = compteur+1;
     end
end
 
 part_coord_x_reshape=zeros(compteur,Var1); % redimensionnement de la matrice contenant les coordonnées x des particules
 part_coord_y_reshape=zeros(compteur,Var1); % redimensionnement de la matrice contenant les coordonnées y des particules
 
 compteur = 0.;
 
 % boucle enregistrant dans les matrices précédentes que les coordonnées des particules visibles de la première à la dernière image
 for i=1:size(part_coord_x,1);
     if (all(part_coord_x(i,:))==1); % ne rentre dans la boucle que si aucune case n'est vide
         compteur = compteur +1;
         part_coord_x_reshape(compteur,:)=part_coord_x(i,:); 
         part_coord_y_reshape(compteur,:)=part_coord_y(i,:);
     end
 end
 
 % représentation graphique des lignes de courant obtenues
 figure(1);
 for i=1:size(part_coord_x_reshape,1);
     plot(part_coord_x_reshape(i,:),part_coord_y_reshape(i,:));
     hold on
 end
 axis image
 saveas(figure(1),VarN,'png') % sauvegarde de la figure
 hold off