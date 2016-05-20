function [varm]=SiStER_interp_shear_nodes_to_markers(varS,x,y,xm,ym,icn,jcn)
% [varm]=SiStER_interp_shear_nodes_to_markers(varS,x,y,xm,ym,icn,jcn)
% interpolates properties from shear nodes to markers

[m, ~] = size(varS);

xnodes = [x(jcn); x(jcn+1); x(jcn+1); x(jcn)];
ynodes = [y(icn); y(icn); y(icn+1); y(icn+1)];

INDEX = sub2ind(size(varS), icn, jcn);
VARnodes = [varS(INDEX); varS(INDEX+m); varS(INDEX+m+1); varS(INDEX+1)];

varm = SiStER_interp_grid_to_marker_vector(xnodes,ynodes,VARnodes,xm,ym);

end

