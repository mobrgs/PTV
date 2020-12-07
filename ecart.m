function [evalecart] = ecart(vect_sol,I,sigma)
%Cette fonction permet de retrouver des coordonn�es en fonction d'une
%r�partition gaussienne des intensit� de niveaux de gris autour de ces
%coordonn�es

[X,Y]=meshgrid(1:size(I,1),1:size(I,2));

evalecart = vect_sol(3)*exp((-((X-vect_sol(1)).^2)+(Y-vect_sol(2)).^2)./2*sigma^2)-I;

end


