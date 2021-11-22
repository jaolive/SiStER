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
    ep1=MAT(types(i)).ep1;    %Start weakening at this plastic strain

    % get cohesion
    %cohes(imInd)=max(Cmax+(Cmin-Cmax).*ep(imInd)./epscrit,Cmin);
    
    % initial weakening cohesion - from G Ito
    % (https://github.com/GTAIto/JAOlive_SiStER) - TMorrow 22 Nov 2021
    cohes(imInd)=min( max( Cmax+(Cmin-Cmax).*(ep(imInd)-ep1)./(epscrit-ep1),Cmin),Cmax);
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


