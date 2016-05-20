% SiStER Initialize


% construct staggered grids
[X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid(xsize,ysize,GRID);

% initialize marker arrays and positions
[xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad);

% locate markers with respect to grid
[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);

% assign marker phases
[im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym);
% initialize marker plastic strain (to zero) and strain rate (to one)
ep=zeros(size(xm));
epNH=ep;
epsIIm=ones(size(xm));

% initialize marker stresses
sxxm=zeros(size(xm));
sxym=sxxm;

% initialize marker index (a unique number to identify and track each marker)
idm=1:length(xm);

% initialize temperature structure on nodes
T=PARAMS.a0+PARAMS.a1*Y+PARAMS.a2*Y.^2+PARAMS.a3*Y.^3;
T=T+PARAMS.amp*sin(2*pi*X/PARAMS.lam);
if PARAMS.ynTreset==1 % reset T=T0 in top layer
    T(T<PARAMS.T0)=PARAMS.T0;
end
% pass initial nodal T to markers
[Tm]=SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
Tm0=Tm;

% initialize nodal strain rate
EXX=zeros(size(X));
EXY=zeros(size(X));
vx=zeros(size(X));
vy=zeros(size(X));

% initialize marker chain to track base of layer 1 (sticky layer)
Ntopo=1000;
topo_x=linspace(0,xsize,Ntopo);
topo_y=GEOM(1).bot*ones(size(topo_x));
topo_marker_spacing=mean(diff(topo_x)); % initial mean spacing of topography markers


