function [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,varargin)
% [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,varargin)
% interpolates marker properties to shear nodes
% First cut - J.A. Olive, March 2011
% Modified by E. Mittelstaedt, April 2011, to allow variable inputs.  
% Modified by B.Z. Klein, Spring 2014, for speedup
% Modified by B.Z. Klein, Summer 2014, for further speedup (vectorized)


Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

% MITTELSTAEDT - check for number of properties to interpolate
numV = size(varargin,2);


% MITTELSTAEDT % establish interpolants matrices
n2interp = repmat(struct('data', zeros(Ny,Nx)), 1, numV);

JCN = interp1(x, 1:length(x), xm, 'nearest', 'extrap'); %% these are like the jcn and icn elsewhere, except the nodes are centered instead of upper left.
ICN = interp1(y, 1:length(y), ym, 'nearest', 'extrap'); %% this makes a lot of the indexing much simpler below.


%% Interior Cells

center = jcn>1 & jcn<Nx & icn>1 & icn<Ny;
shiftLeft = jcn<Nx-1 & icn>1 & icn<Ny;
shiftUp = jcn>1 & jcn<Nx & icn<Ny-1;
shiftBoth = jcn<Nx-1 & icn<Ny-1;


cell1 = center & ((xm-x(JCN)) > 0) & ((ym - y(ICN)) > 0);  %% these are logical arrays that index the original quadrants
cell2 = shiftLeft & ((xm-x(JCN)) < 0) & ((ym - y(ICN)) > 0);
cell3 = shiftBoth & ((xm-x(JCN)) < 0) & ((ym - y(ICN)) < 0);
cell4 = shiftUp & ((xm-x(JCN)) > 0) & ((ym - y(ICN)) < 0);

%%% WEIGHTING (equal for now because that is what I'm running)

wc1 = 0.25;
wc2 = 0.25;
wc3 = 0.25;
wc4 = 0.25;


% cell 1 (i,j,1)

dxm = xm(cell1) - x(JCN(cell1));
dym = ym(cell1) - y(ICN(cell1));
ddx = dx(JCN(cell1));
ddy = dy(ICN(cell1));

wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', JCN(cell1)'], wm1);

% cell 2 (i, j-1, 2)

dxm = xm(cell2) - x(JCN(cell2));
dym = ym(cell2) - y(ICN(cell2));
ddx = dx(JCN(cell2)-1);
ddy = dy(ICN(cell2));

wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2 = accumarray([ICN(cell2)', JCN(cell2)'], wm2);

% cell 3 (i-1, j-1, 3)

dxm = xm(cell3) - x(JCN(cell3));
dym = ym(cell3) - y(ICN(cell3));
ddx = dx(JCN(cell3)-1);
ddy = dy(ICN(cell3)-1);

wm3 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w3 = accumarray([ICN(cell3)', JCN(cell3)'], wm3);

% cell 4 (i-1, j, 4)

dxm = xm(cell4) - x(JCN(cell4));
dym = ym(cell4) - y(ICN(cell4));
ddx = dx(JCN(cell4));
ddy = dy(ICN(cell4)-1);

wm4 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w4 = accumarray([ICN(cell4)', JCN(cell4)'], wm4);

%loop over material properties to interpolate

for vn = 1:numV
    n2interp(vn).data = (wc1*accumarray([ICN(cell1)', JCN(cell1)'], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', JCN(cell2)'], varargin{vn}(cell2).*wm2)./w2 + ...
        wc3*accumarray([ICN(cell3)', JCN(cell3)'], varargin{vn}(cell3).*wm3)./w3 + ...
        wc4*accumarray([ICN(cell4)', JCN(cell4)'], varargin{vn}(cell4).*wm4)./w4)./...
        (wc1+wc2+wc4+wc4);
end



%% EDGES

%%% top edge

topEdge = jcn>1 & jcn<Nx & icn==1;
shifted = jcn<Nx-1 & icn==1;

% cell 1

cell1 = shifted & quad==2;

ddx = dx(JCN(cell1)-1);
ddy = dy(1);
dxm = xm(cell1) - x(JCN(cell1));
dym = ym(cell1) - y(ICN(cell1));
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', JCN(cell1)'], wm1);

% cell 2 

cell2 = topEdge & quad==1;

ddx = dx(JCN(cell2));
ddy = dy(1);
dxm = xm(cell2) - x(JCN(cell2));
dym = ym(cell2) - y(ICN(cell2));
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2  = accumarray([ICN(cell2)', JCN(cell2)'], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ICN(cell1)', JCN(cell1)'], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', JCN(cell2)'], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(1,2:end) = temp(2:end);
end

clear w1 w2

%%% bottom edge

bottomEdge = jcn>1 & jcn<Nx & icn==Ny-1;
shifted =    jcn<Nx-1       & icn==Ny-1;

% cell 1

cell1 = shifted & quad==3;

ddx = dx(JCN(cell1)-1);
ddy = dy(Ny-1);
dxm = xm(cell1) - x(JCN(cell1));
dym = ym(cell1) - y(end-1);
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ones(sum(cell1),1), JCN(cell1)'], wm1);

% cell 2

cell2 = bottomEdge & quad==4;

ddx = dx(JCN(cell2));
ddy = dy(Ny-1);
dxm = xm(cell2) - x(JCN(cell2));
dym = ym(cell2) - y(end-1);
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2  = accumarray([ones(sum(cell2),1), JCN(cell2)'], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ones(sum(cell1),1), JCN(cell1)'], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ones(sum(cell2),1), JCN(cell2)'], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(Ny,2:end) = temp(2:end);
end

%%% left edge

leftEdge = jcn==1 & icn>1 & icn<Ny;
shifted  = jcn==1 & icn<Ny-1;

% cell 1

cell1 = shifted & quad==4;

ddx = dx(1);
ddy = dy(ICN(cell1)-1);
dxm = xm(cell1) - x(1);
dym = ym(cell1) - y(ICN(cell1));
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', ones(sum(cell1),1)], wm1);

% cell 2

cell2 = leftEdge & quad==1;

ddx = dx(1);
ddy = dy(ICN(cell2));
dxm = xm(cell2) - x(1);
dym = ym(cell2) - y(ICN(cell2));
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ICN(cell1)', ones(sum(cell1),1)], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', ones(sum(cell2),1)], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(2:end-1, 1) = temp(2:end);
end

%%% right edge

rightEdge = jcn==Nx-1 & icn>1 & icn<Ny;
shifted =   jcn==Nx-1 & icn<Ny-1;

% cell 1

cell1 = shifted & quad==3;

ddx = dx(Nx-1);
ddy = dy(ICN(cell1)-1);
dxm = xm(cell1) - x(Nx-1);
dym = ym(cell1) - y(ICN(cell1));
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', ones(sum(cell1),1)], wm1);

% cell 2

cell2 = rightEdge & quad==2;

ddx = dx(Nx-1);
ddy = dy(ICN(cell2));
dxm = xm(cell2) - x(Nx-1);
dym = ym(cell2) - y(ICN(cell2));
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ICN(cell1)', ones(sum(cell1),1)], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', ones(sum(cell2),1)], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(2:end-1, Nx) = temp(2:end);
end

%% CORNERS

% upper left

upperLeft = jcn==1 & icn==1 & quad==1;

ddx = dx(1);
ddy = dy(1);
dxm = xm(upperLeft) - x(1);
dym = ym(upperLeft) - y(1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(1,1) = sum(varargin{vn}(upperLeft).*wm)./wco;
end

% upper right

upperRight = icn==1 & jcn==Nx-1 & quad==2;

ddx = dx(Nx-1);
ddy = dy(1);
dxm = xm(upperRight) - x(Nx-1);
dym = ym(upperRight) - y(1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(1,Nx) = sum(varargin{vn}(upperRight).*wm)./wco;
end

% lower Right

lowerRight = icn==Ny-1 & jcn==Nx-1 & quad==3;

ddx = dx(Nx-1);
ddy = dy(Ny-1);
dxm = xm(lowerRight) - x(Nx-1);
dym = ym(lowerRight) - y(Ny-1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(Ny,Nx) = sum(varargin{vn}(lowerRight).*wm)./wco;
end

% lower left

lowerLeft = icn==Ny-1 & jcn==1 & quad==4;

ddx = dx(1);
ddy = dy(Ny-1);
dxm = xm(lowerLeft) - x(1);
dym = ym(lowerLeft) - y(Ny-1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(Ny,1) = sum(varargin{vn}(lowerLeft).*wm)./wco;
end


