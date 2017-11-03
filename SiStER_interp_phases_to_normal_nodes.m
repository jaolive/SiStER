function [phaseWeights] = SiStER_interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases, maxPhases)
% [phaseWeights] = SiStER_interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases, maxPhases)
% B.Z. Klein, July 2017, an interp function specific to phase (for normal nodes), to enable 
% exact mixing of several phases

Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);


INDEX = sub2ind([Ny-1, Nx-1], icn, jcn);

AcCell = bsxfun(@times, dy', dx);

xN = x(1:Nx-1) + dx/2;
yN = y(1:Ny-1) + dy/2;
[XN, YN] = meshgrid(xN, yN);

AMvec = abs((xm - XN(INDEX)).*(ym - YN(INDEX)));
WMvec = (AcCell(INDEX) - AMvec)./AcCell(INDEX);

phaseWeights=zeros(Ny,Nx,maxPhases);

for n = 1:maxPhases
    
    phaseMask = phases==n;
    phaseWeights(2:Ny,2:Nx,n) = accumarray([icn(phaseMask)' jcn(phaseMask)'], WMvec(phaseMask)', size(AcCell), @sum);
    
end
   
sumWeights = repmat(sum(phaseWeights, 3), [1, 1, maxPhases]);
phaseWeights = phaseWeights./sumWeights;

end




