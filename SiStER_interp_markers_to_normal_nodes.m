function [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)
% [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)
% interpolates properties (in the order of input) from markers to normal nodes
%
% First cut - J.A. Olive, March 2011
% Modified by E. Mittelstaedt, April 2011 to allow multiple inputs  
% Modified by B.Z. Klein, Spring 2014 for speedup

Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

% check for number of properties to interpolate
numV = size(varargin,2);
n2interp(1:numV) = struct('data',zeros(Ny,Nx));

INDEX = sub2ind([Ny-1, Nx-1], icn, jcn);

AcCell = bsxfun(@times, dy', dx);

xN = x(1:Nx-1) + dx/2;
yN = y(1:Ny-1) + dy/2;
[XN, YN] = meshgrid(xN, yN);

AMvec = abs((xm - XN(INDEX)).*(ym - YN(INDEX)));
WMvec = (AcCell(INDEX) - AMvec)./AcCell(INDEX);

w = accumarray([icn' jcn'], WMvec', [], @sum);

for vn = 1:numV
    VecData = varargin{vn}.*WMvec;
    n2interp(vn).data(2:Ny,2:Nx) = accumarray([icn' jcn'], VecData', [], @sum)./w;
end

end




