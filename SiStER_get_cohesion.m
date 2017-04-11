function [cohes]=SiStER_get_cohesion(im,ep,MAT)
% [cohes]=SiStER_get_cohesion(im,ep,MAT)
% compute cohesion on markers based on ep
%G.Ito 8/2016
% sped up by B. Klein 9/2016

cohes = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    Cmax=MAT(types(i)).Cmax;
    Cmin=MAT(types(i)).Cmin;
    epscrit=MAT(types(i)).ecrit;

    % get cohesion
    cohes(imInd)=max(Cmax+(Cmin-Cmax).*ep(imInd)./epscrit,Cmin);
end

return

%% OLD VERSION the bracket approach is really slow
% 
% Cmax=[MAT(im).Cmax];
% Cmin=[MAT(im).Cmin];
% epscrit=[MAT(im).ecrit];
% 
% % get cohesion
% cohes=max(Cmax+(Cmin-Cmax).*ep./epscrit,Cmin);
% 
% return


