%=========================================================================
function [varS]=SiStER_interp_normal_to_shear_nodes(varN,dx,dy)
% Interpolates values on normal nodes (varN) to shear nodes
% Potential shortcoming:  interpolation to side shear nodes are only
% based on the closest normal nodes and thus values will be horizontally
% uniform with 1/2 cell of the left/right sides and vertically uniform
% with 1/2 cell of the top and bottoms
%
% G.Ito 8/2016
%=========================================================================
varS=zeros(size(varN));
[Ny,Nx]=size(varS);
i1=2:Ny-1; j1=2:Nx-1; 
i2=3:Ny;   j2=3:Nx;
dydx=varS;
dydx(2:end,2:end)=dy'*dx;

%--------------------------------------------------------------------------
% interior nodes
%--------------------------------------------------------------------------
varS(i1,j1)=(varN(i1,j1).*dydx(i2,j2)/4 + varN(i1,j2).*dydx(i2,j1)/4 +...
             varN(i2,j1).*dydx(i1,j2)/4 + varN(i2,j2).*dydx(i1,j1)/4)./...
             ((dydx(i1,j1)+dydx(i1,j2)+dydx(i2,j1)+dydx(i2,j2))./4);
                   
%--------------------------------------------------------------------------
% top and bottom, excluding corners
%--------------------------------------------------------------------------
varS(1,j1)= (varN(2,j1).*dydx(2,j2)/4 + varN(2,j2).*dydx(2,j1)/4)./...
            (dydx(2,j2)/4 + dydx(2,j1)/4);
                       
varS(Ny,j1)=(varN(Ny,j1).*dydx(Ny,j2)/4 + varN(Ny,j2).*dydx(Ny,j1)/4)./...
            (dydx(Ny,j2)/4 + dydx(Ny,j1)/4);
                        
%--------------------------------------------------------------------------
% left and right, excluding corners
%--------------------------------------------------------------------------
varS(i1,1)= (varN(i1,2).*dydx(i2,2)/4 + varN(i2,2).*dydx(i1,2)/4)./...
            (dydx(i1,2)/4 + dydx(i2,2)/4);
        
varS(i1,Nx)=(varN(i1,Nx).*dydx(i2,Nx)/4 + varN(i2,Nx).*dydx(i1,Nx)/4)./...
            (dydx(i2,Nx)./4 +dydx(i1,Nx)/4);
        
%--------------------------------------------------------------------------
% corners
%--------------------------------------------------------------------------
varS(1,1) = varN(2,2); varS(1,Nx) =varN(2,Nx);
varS(Ny,1)=varN(Ny,2); varS(Ny,Nx)=varN(Ny,Nx);
                   
                   
                      
return
