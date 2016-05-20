function [rhom]=SiStER_get_density(im,Tm,MAT)
% [rho]=SiStER_get_density(im,Tm,MAT)
% obtain density from temperature and material identity
% 

T0=0;
rhom = [MAT(im).rho0].*(1-[MAT(im).alpha].*(Tm-T0));