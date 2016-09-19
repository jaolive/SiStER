function [cohes]=SiStER_get_cohesion(im,ep,MAT)
% [cohes]=SiStER_get_cohesion(im,ep,MAT)
% compute cohesion on markers based on ep
%G.Ito 8/2016

Cmax=[MAT(im).Cmax];
Cmin=[MAT(im).Cmin];
epscrit=[MAT(im).ecrit];

% get cohesion
cohes=max(Cmax+(Cmin-Cmax).*ep./epscrit,Cmin);

return