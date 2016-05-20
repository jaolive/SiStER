function [ys]=SiStER_get_yield_stress(im,ep,p,MAT)
% [ys]=SiStER_get_yield_stress(im,ep,p,MAT)
% obtain yield stress from material identity, pressure 
% and accumulated plastic strain

% Mohr-Coulomb written as Drucker-Prager:
% sII = A*P+B;
% where A = sind(atand(mu))
% and B = C*cosd(atand(mu))

mu=[MAT(im).mu];
Cmax=[MAT(im).Cmax];
Cmin=[MAT(im).Cmin];
epscrit=[MAT(im).ecrit];

% get cohesion
cohes=max(Cmax+(Cmin-Cmax).*ep./epscrit,Cmin);
% cohesion weakens linearly with accumulated plastic strain
% until minimum cohesion is reached (e.g., Lavier et al., 2000 JGR)

ys=cohes.*cosd(atand(mu))+sind(atand(mu)).*p; % Mohr-Coulomb

