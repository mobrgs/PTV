function [evalecart] = ecart(vect_sol,I,sigma)
%Cette fonction permet de retrouver des coordonnées en fonction d'une
%répartition gaussienne des intensité de niveaux de gris autour de ces
%coordonnées

[X,Y]=meshgrid(1:size(I,1),1:size(I,2));

evalecart = vect_sol(3)*exp((-((X-vect_sol(1)).^2)+(Y-vect_sol(2)).^2)./2*sigma^2)-I;

end


