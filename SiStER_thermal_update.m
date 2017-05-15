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



% GET VARIABLE DIFFUSIVITY AND CP
if isfield(MAT,'cp')==0 || isfield(MAT,'k')==0   
    cpfield=PARAMS.cpref*ones(size(T));
    kfield=PARAMS.kref*ones(size(T));
    rhofield=PARAMS.rhoref*ones(size(T));    
else
    [km, cpm]=SiStER_get_thermal_properties(im,MAT);
    [rhom]=SiStER_get_density(im,Tm,MAT);
    [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,km,cpm,rhom);
    kfield=n2interp(1).data;
    cpfield=n2interp(2).data;
    rhofield=n2interp(3).data;

end
% THERMAL SOLVE
[T]=SiStER_thermal_solver_sparse_CFD(x,y,Told,rhofield,cpfield,kfield,dt_m,BCtherm,zeros(size(T)));



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


