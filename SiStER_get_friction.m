function [fric]=SiStER_get_friction(im,ep,MAT)
% [fric]=SiStER_get_friction(im,ep,MAT)
% compute friction on markers based on plastic strain, like cohesion
% J.-A. Olive 4/2017


fric = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    mumax=MAT(types(i)).mu;

    if isfield(MAT,'mumin')==1  
        mumin=MAT(types(i)).mumin;
    else
        mumin=MAT(types(i)).mu;
    end

    epscrit=MAT(types(i)).ecrit;

    % get cohesion
    fric(imInd)=max(mumax+(mumin-mumax).*ep(imInd)./epscrit,mumin);
end

return


%% OLD (slow due to brackets)
% mumax=[MAT(im).mu];
% 
% 
% if isfield(MAT,'mumin')==1  
%     mumin=[MAT(im).mumin];
% else
%     mumin=[MAT(im).mu];
% end
% 
% epscrit=[MAT(im).ecrit];
% 
% % get cohesion
% fric=max(mumax+(mumin-mumax).*ep./epscrit,mumin);
% 
% return