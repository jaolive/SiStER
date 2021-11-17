
%=========================================================================
% PLASTIC (EXCLUDING ELASTIC) STRAIN ACCUMULATION
% Not sure whether to use the actual non-elastic strain (i.e., with
% current stresses) or the approximate, in which stress is assumed to equal
% yield stress
% G.Ito 8/16; JAO 9/15 for non-healed ep, fixed by JAO 4/2017
%=========================================================================

dep_s=zeros(size(epsII_s));

if (BC.DIKE.on==1) % epsII_s corrected by removing dike strain when dike is activated
                    % after G.Ito (https://github.com/GTAIto/JAOlive_SiStER) T.Morrow 16 Nov 2021
    EXX_s=SiStER_interp_normal_to_shear_nodes(EXX-BC.DIKE.MV./BC.DIKE.DX,dx,dy);
    dum=sqrt(EXX_s.^2+EXY.^2);  %This is epsII_s
    dep_s(s_nodes_yield) = dt_m.*max(dum(s_nodes_yield)-...
(yield_s(s_nodes_yield)-sqrt(sxxOLD_s(s_nodes_yield).^2+sxyOLD(s_nodes_yield).^2))./(2.*Gs(s_nodes_yield).*dt_m),min(dum(:))*1e-6);
else
    dep_s(s_nodes_yield) = dt_m.*max(epsII_s(s_nodes_yield)-...
(yield_s(s_nodes_yield)-sqrt(sxxOLD_s(s_nodes_yield).^2+sxyOLD(s_nodes_yield).^2))./(2.*Gs(s_nodes_yield).*dt_m),min(epsII_s(:))*1e-6);
end

[depm]=SiStER_interp_shear_nodes_to_markers(reshape(dep_s,Ny,Nx),x,y,xm,ym,icn,jcn);

ep=(ep+depm)./(dt_m/PARAMS.tau_heal+1);
epNH=epNH+depm;