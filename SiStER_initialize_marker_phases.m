function [im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym)
% [im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym)
% this is where the identity of each marker (e.g., air, crust, mantle...)
% gets assigned following the geometry specified in the input file 
% (flat layers, circle or rectangle)
% an alternative geometry can be added here.

% assign material identity on markers
im=zeros(size(xm));


for kk = 1:Nphase
    
    if GEOM(kk).type==1 % layer
        
        im(ym>=GEOM(kk).top & ym<GEOM(kk).bot)=kk;
        
    elseif GEOM(kk).type==2 % circular inclusion
        
        rm=sqrt((xm-GEOM(kk).x0).^2 + (ym-GEOM(kk).y0).^2);
        im(rm<=GEOM(kk).rad)=kk;
        
    elseif GEOM(kk).type==3 % rectangle
        
        im(ym>=GEOM(kk).top & ym<GEOM(kk).bot & xm>=GEOM(kk).left & xm<GEOM(kk).right)=kk;
        
    end
    
end
        
