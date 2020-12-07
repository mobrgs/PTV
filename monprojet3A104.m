% mon projet Matlab de detection / suivi de particules
clear % efface le contenu du Workspace
close all % ferme toutes les figures
clc % efface le contenu de la fenetre de travail (Command Window)
%==========================================================================
% Parametres utilisateur
Var1 = 100; % Nombre d'image � traiter
VarN = 'NomFigure'; % Titre de la figure des r�sultats enregistr�e sous forme .png
% Nothing to modify below this line
%==========================================================================
% filename : structure contenant les noms des images a traiter
% pathname : chemin d'acces au repertoire contenant les images a traiter
% C : vecteur contenant les coordonn�es colonne des pixels des maximums locaux (image amont)
% L : vecteur contenant les coordonn�es ligne des pixels des maximums locaux (image amont)
% C2 : vecteur contenant les coordonn�es colonne des pixels des maximums locaux (image avale)
% L2 : vecteur contenant les coordonn�e ligne des pixels des maximums locaux (image avales)
% C2_zoom : vecteur contenant les coordonn�es colonne des pixels des maximums locaux dans une partie zoom�e de l'image avale
% C2_zoom : vecteur contenant les coordonn�es ligne des pixels des maximums locaux dans une partie zoom�e de l'image avale
% dilate_I : image amont dont les maximums locaux sont �tendus
% dilate_I2 : image avale dont les maximums locaux sont �tendus
% dilate_zoom_I2 : partie zoom�e de l'image avale dont les maximums locaux sont �tendus
% ext_zoom_I2 : zoom_I2 compl�t� de pixels noirs sur se p�riph�rie de taille "len"
% I1 : matrice contenant les intensit�s des pixels de l'image amont
% I2 : matrice contenant les intensit� des pixels de l'image avale
% I1_ext : I1 �tendue d'une taille "len" de lignes et colonnes nulles
% I2_ext : I2 �tendue d'une taille "len" de lignes et colonnes nulles
% mat_corr : matrice contenant les valeurs de corr�lation entre les particules de l'image amont (lignes) et avales (colonnes)
% mat_sol : matrice dont la premi�re colonne prend le num�ro de particule de l'image amont et la deuxi�me colonne les num�ro de particule de l'image avale correspondante
% N : taille du zoom effectiu� de part et d'autre de chaque pixel pris comme maximum local
% N2 : taille du zoom dans la zone de d�placement d'une particule d'une image � l'autre
% part_coord_x : coordonn�es x des centres des particules (lignes : num�ro particule ; colonne : num�ro image)
% part_coord_y : coordonn�es y des centres des particules (lignes : num�ro particule ; colonne : num�ro image)
% part_init : particule regard�e pour la corr�lation sur l'image amont
% part_test : particule regard�e pour la corr�lation sur l'image avale
% pos : vecteur contenant les position lin�aire des maximums locaux d'une image
% R : coefficient de corr�lation
% seuil : seuil sous lequel les maximums locaux ne sont pas d�tect�s (image enti�re)
% seuil2 : seuil sous lequel les maximums locaux ne sont pas d�tect�s (partie zoom�e de l'image avale)
% sigma : �cart type de la fonction gaussienne
% Xsol : coordonn�es x des positions exactes des paticules sur l'image amont
% Ysol : coordonn�es y des positions exactes des paticules sur l'image amont
% Xsol2 : coordonn�es x des positions exactes des paticules sur l'image avale
% Ysol2 : coordonn�es y des positions exactes des paticules sur l'image avale
% zoom_I2 : partie zoom�e de l'image avale pour la r�cherche des particule correspondantes � celles pr�sente dans l'image amont
% zoom : zoom autour d'une particule � la distance N du pixel central
% compteur : pour contenant le nombre de fois o� une boucle est ex�cut�e

[filename, pathname] = uigetfile({' * .tiff; * .png'}, 'Select ...images','MultiSelect','on');
if isequal(filename,0)
disp('User selected Cancel')
return
end
cd(pathname); % rend le r�pertoire d'acces des images actif

load('parameters.mat') % charge les param�tres du probl�me

%% Lecture de la premi�re image

 sigma = 1.;
 I = imread(filename{1}); % conversion de l'image en matrice remplie de son intensit� en niveau de gris par pixel
 
%% Boucle d'initialisation pour avoir les coordonn�es des particules � traquer depuis l'image 1
 
     I = imread(filename{1}); % conversion de l'image 1 en matrice remplie de son intensit� en niveau de gris
     I=im2double(I); % mise en double pr�cision des valeurs contenues dans la matrice
     len = 5; 
     I_ext = wextend(2,'zpd',I,len); % �tendue de la matrice avec des ligne et colonne nulles sur une taille "len"
     
     SE = strel('disk', 3);
     dilate_I=imdilate(I_ext,SE); % �tendue des intensit�s des maximums locaux sur la taille estim�e d'une particule

     seuil = multithresh(dilate_I); % estimation du seuil pour rep�rer les maximums locaux
     pos=find((I_ext)==(dilate_I) & (I_ext)>=seuil); % d�termination des positions lin�aire de ces maximums locaux
     [L,C]=ind2sub(size(I_ext), pos); % convertion de ces positions lin�aires en coordonn�es ligne, colonne

     part_coord_x = zeros(size(L,1),Var1); % cr�ation de la matrice contenant les coordonn�es x des particules pour toutes les images
     part_coord_y = zeros(size(L,1),Var1); % cr�ation de la matrice contenant les coordonn�es y des particules pour toutes les images

     N=2; 

     Xsol=zeros(size(C,1),1); % cr�ation du vecteur solution des coordonn�es x des particules pour l'image 1
     Ysol=zeros(size(L,1),1); % cr�ation du vecteur solution des coordonn�es x des particules pour l'image 1

     % boucle permettant de calculer la position exacte des centres des particules
     for i = 1:size(L,1);
        zoom=I_ext(L(i)-N:L(i)+N,C(i)-N:C(i)+N);
        vect_sol=[0 0 0];
        [Imax]=max(zoom(:));
        [solution]=lsqnonlin(@(vect_sol)ecart(vect_sol,zoom,sigma),[N, N, Imax]); % fonction de minimisation en concid�rant une r�partition gaussienne de l'intensit� autour du centre
         Xsol(i)=solution(1)+C(i)-N-1;
         Ysol(i)=solution(2)+L(i)-N-1;
     end
     
     % boucle remplissant la premi�re colonne de la matrice enregistrant les coordonn�es des particules
     % ce vont �tre les particules traqu�es
     for i=1:size(part_coord_x,1);
        part_coord_x(i,1)=Xsol(i);
        part_coord_y(i,1)=Ysol(i);
     end

%% Boucle sur toutes les images

for n=1:Var1-1;

    % boucle transferant les coordonn�es des maximums locaux pr�c�demment enregistr�es pour l'image avale comme les coordonn�es de l'images amont
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
    
     I1 = imread(filename{n}); % conversion de l'image amont en matrice remplie de son intensit� en niveau de gris
     I2 = imread(filename{n+1}); % conversion de l'image avale en matrice remplie de son intensit� en niveau de gris

     I1=im2double(I1); % mise en double pr�cision des valeurs contenues dans la matrice
     I2=im2double(I2); % mise en double pr�cision des valeurs contenues dans la matrice

     I1_ext = wextend(2,'zpd',I1,len); % �tendue de la matrice avec des ligne et colonne nulles sur une taille "len"
     I2_ext = wextend(2,'zpd',I2,len); % �tendue de la matrice avec des ligne et colonne nulles sur une taille "len"

     SE = strel('disk', 3);
     dilate_I2=imdilate(I2_ext,SE); % �tendue des intensit�s des maximums locaux sur la taille estim�e d'une particule de l'image avale
     seuil = multithresh(dilate_I2); % estimation du seuil pour rep�rer les maximums locaux
     pos=find((I2_ext)==(dilate_I2) & (I2_ext)>=seuil); % d�termination des positions lin�aire de ces maximums locaux
     [L2,C2]=ind2sub(size(I2_ext), pos); % convertion de ces positions lin�aires en coordonn�es ligne, colonne

     mat_corr=zeros(size(L,1), size(L2,1)); % cr�ation de la matrice de corr�lation

     N2=5;

    % boucle calculant la corr�lation entre une particule de l'image amont et des particules de l'image avale contenunes dans une zone pr�d�finies
    for j=1:size(L,1);
        if (L(j)~=0);
         zoom_I2=I2_ext(L(j)-N2:L(j)+N2,C(j)-N2:C(j)+N2); % d�finition de la zone de recherche sur l'image avale
         SE2 = strel('disk', 3);
         % recherche des maximums locaux dans cette zone pr�d�finies
         dilate_zoom_I2=imdilate(zoom_I2,SE2);  % �tendue des intensit�s des maximums locaux sur la taille estim�e d'une particule
         seuil2 = multithresh(dilate_zoom_I2);  % estimation du seuil pour rep�rer les maximums locaux
         pos2=find((zoom_I2)==(dilate_zoom_I2) & (zoom_I2)>=seuil2);  % d�termination des positions lin�aire de ces maximums locaux
         [L2_zoom,C2_zoom]=ind2sub(size(zoom_I2), pos2); % convertion de ces positions lin�aires en coordonn�es ligne, colonne

         part_init=I1_ext(L(j)-N:L(j)+N,C(j)-N:C(j)+N); % d�coupage autour de la particule de l'image amont
         ext_zoom_I2=wextend(2,'zpd',zoom_I2,len);

             % calcul de la coor�lation de toutes les particules pr�sentes dans la zone pr�d�finie pr�c�demment avec la particule de l'image amont
             for k=1:size(C2_zoom,1)
                 part_test=ext_zoom_I2(L2_zoom(k)+len-N:L2_zoom(k)+len+N,C2_zoom(k)+len-N:C2_zoom(k)+len+N); % d�coupage autour de la particule � tester dans la zone pr�d�finie de l'image avale
                 R=corr2(part_init,part_test); % calcul du coefficient de cor��lation entre les deux particules
                 num=find((C2==C2_zoom(k)+C(j)-N2-1)&(L2==L2_zoom(k)+L(j)-N2-1)); % recherche du num�ro de la particule de l'image avale
                 mat_corr(j,num)=R; % remplissage de la matrice de corr�lation (ligne : num�ro particule image amont, colonne : num�ro particule image avale)
             end
        end
    end

    % cr�ation de la matrice de correspondance entre les particules de l'image amont et de l'image avale
    mat_sol=zeros(min(size(mat_corr)),2);

    % boucle remplissant cette matrice de correspondance
    for k=1:size(mat_sol,1);
        if (any(mat_corr(k,:))==1);
            a=find(mat_corr(k,:)==max(mat_corr(k,:))); % recherche de la corr�lation maximum entre les particules des deux images
            if ((mat_corr(k,a))==max(mat_corr(:,a))); % v�rification de l'unicit� de la solution
                mat_sol(k,1)=k; % remplissage de la premi�re colonne de la matrice de correspondance
                mat_sol(k,2)=a; % remplissage de la deuxi�me colonne de la matrice de correspondance
            end
        end
    end

     Xsol2=zeros(size(C2,1),1); % cr�ation du vecteur enregistrant les coordonn�es x des centres des particules
     Ysol2=zeros(size(L2,1),1); % cr�ation du vecteur enregistrant les coordonn�es y des centres des particules

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
            part_coord_x(i,n+1)=Xsol2(mat_sol(i,2)); % coordonn�es x
            part_coord_y(i,n+1)=Ysol2(mat_sol(i,2)); % coordonn�es y
        end
    end 
 end

 % save("Projet_donnees"); % sauvegarde le workspace

%% Ecriture du graphique pr�sentant les r�sultats de la m�thode PTV

compteur = 0.; 

% boucle comptant le nombre de particule visible sur toutes les images
for i=1:size(part_coord_x,1);
     if (all(part_coord_x(i,:))==1); % ne rentre dans la boucle que si aucune case n'est vide
         compteur = compteur+1;
     end
end
 
 part_coord_x_reshape=zeros(compteur,Var1); % redimensionnement de la matrice contenant les coordonn�es x des particules
 part_coord_y_reshape=zeros(compteur,Var1); % redimensionnement de la matrice contenant les coordonn�es y des particules
 
 compteur = 0.;
 
 % boucle enregistrant dans les matrices pr�c�dentes que les coordonn�es des particules visibles de la premi�re � la derni�re image
 for i=1:size(part_coord_x,1);
     if (all(part_coord_x(i,:))==1); % ne rentre dans la boucle que si aucune case n'est vide
         compteur = compteur +1;
         part_coord_x_reshape(compteur,:)=part_coord_x(i,:); 
         part_coord_y_reshape(compteur,:)=part_coord_y(i,:);
     end
 end
 
 % repr�sentation graphique des lignes de courant obtenues
 figure(1);
 for i=1:size(part_coord_x_reshape,1);
     plot(part_coord_x_reshape(i,:),part_coord_y_reshape(i,:));
     hold on
 end
 axis image
 saveas(figure(1),VarN,'png') % sauvegarde de la figure
 hold off