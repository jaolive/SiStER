function [Gm]=SiStER_get_elastic_moduli(im,MAT)
% [ys]=SiStER_get_elastic_moduli(im,MAT)
% obtain shear modulus from material properties
% B. Klein 9/16

Gm = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    Gm(imInd) = MAT(types(i)).G;
end


% faster than:
% Gm=[MAT(im).G];



