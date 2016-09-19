%==========================================================================
function [var] = SiStER_get_mean_material_prop(matp,phase_node)
%
% Computes mean value of MAT.{property} at nodes
% e.g., [Mu_s]=SiStER_material_prop(cell2mat({MAT([1:Nphase]).mu}),phase_s);
%
% Shortcoming:  phase_node (=phase_s or phase_n) is phase number determined
% by interpolating im from marker to (shear or normal) node.  Thus this will
% work only if there are two phases adjacent to the nodes and those two 
% phases are consecutively numbered. 
% G.Ito 8/16
%==========================================================================


phase1=floor(phase_node);
var=matp(phase1);

ii=find(phase_node > phase1);
if (length(ii)>1);
    if (sum(ceil(phase_node)>phase1+1));
        disp('problems in SiStER_material_prop_s');
        halt
    end;
   
    f1=((phase1(ii)+1)-phase_node(ii))';  %fraction of phase 1
    var(ii)=f1.*matp(phase1(ii))+(1-f1).*matp(phase1(ii)+1);
end

