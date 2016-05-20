% SiStER THERMAL SOLVE

% get previous temperature on nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
Told(:,:)  = n2interp(1).data;


% enforce Dirichlet boundary conditions to avoid mismatch between markers
% and nodes
if BCtherm.top(1)==1
    Told(1,:)=BCtherm.top(2);
end
if BCtherm.bot(1)==1
    Told(Ny,:)=BCtherm.bot(2);
end
if BCtherm.left(1)==1
    Told(:,1)=BCtherm.left(2);
end
if BCtherm.right(1)==1
    Told(:,Nx)=BCtherm.right(2);
end


% In a future release we will have spatially variable diffusivity
% using conductivity, density and cp carried on markers
% for now we use a constant diffusivity determined from the reference
% values in the input file


[T]=SiStER_thermal_solver_sparse(x,y,Told,PARAMS.rhoref*ones(size(T)),PARAMS.cpref*ones(size(T)),PARAMS.kref*ones(size(T)),PARAMS.kref*ones(size(T)),dt_m,BCtherm);  


% temperature change
dT=T-Told;
% enforce Dirichlet boundary conditions to avoid mismatch between markers
% and nodes
if BCtherm.top(1)==1
    dT(1,:)=0;
end
if BCtherm.bot(1)==1
    dT(Ny,:)=0;
end
if BCtherm.left(1)==1
    dT(:,1)=0;
end
if BCtherm.right(1)==1
    dT(:,Nx)=0;
end

[Tm]=SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);

if PARAMS.ynTreset==1 % reset T=T0 in top layer
     Tm(im==1)=PARAMS.T0;
end


