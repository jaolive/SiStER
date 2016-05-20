function [param] = SiStER_interp_grid_to_marker_vector(xnodes,ynodes,paramnodes,xm,ym)
% [param] = SiStER_interp_grid_to_marker_vector(xnodes,ynodes,paramnodes,xm,ym)
% interpolates a grid value to the entire vector of markers,
% used by SiStER_get_marker_velocities
% B.Z. Klein, 2013

dxm=(xm-xnodes(1,:));
dym=(ym-ynodes(1,:));
dx=abs(xnodes(2,:)-xnodes(1,:));
dy=abs(ynodes(3,:)-ynodes(2,:));

w1=(1-dxm./dx).*(1-dym./dy);
w2=(dxm./dx).*(1-dym./dy);
w4=(dym./dy).*(1-dxm./dx);
w3=(dxm./dx).*(dym./dy);
wnodes=[w1; w2; w3; w4];
    
param=sum(paramnodes.*wnodes, 1);

end