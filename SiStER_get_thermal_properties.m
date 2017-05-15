function [km, cpm]=SiStER_get_thermal_properties(im,MAT)
% [km, cpm]=SiStER_get_thermal_properties(im,MAT)
% obtain thermal conductivity and heat capacity

km = zeros(size(im));
cpm = km;

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    km(imInd) = MAT(types(i)).k;
    cpm(imInd) = MAT(types(i)).cp;
end

