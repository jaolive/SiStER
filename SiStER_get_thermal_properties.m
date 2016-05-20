function [km, cpm]=SiStER_get_thermal_properties(im,MAT)
% [ys]=SiStER_get_thermal_properties(im,MAT)
% obtain thermal conductivity and heat capacity

km=[MAT(im).k];
cpm=[MAT(im).cp];


