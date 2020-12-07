function [niveaugris] = gauss2D(x,y,x0,y0,sigma,Amax)
%Fonction de Gauss 2D qui retourne une intensité en niveau de gris
%   x0 et y0 sont les coordonnées du centre de la particule ; sigma l'écart
%   type, Amax le niveau de gris max (256 en 8 bits)
niveaugris = Amax*exp(-(((x-x0).^2)+((y-y0).^2))./2*sigma^2);

end

