% SiStER_Initialize_dike

% Initializes required dike arrays and property markers
% TMorrow 23 Sep 2019

BC.DIKE.M=zeros(size(Y));
BC.DIKE.MV=zeros(size(Y));
BC.DIKE.DX=repmat([1 dx],numel(dy)+1,1); % j-1 node is the normal-node associated with M
BC.DIKE.DY=repmat([dy 1]',1,numel(dx)+1);
BC.DIKE.H=zeros(size(Y));

% ... and more as needed

% GET MARKER DENSITIES 
[rhom]=SiStER_get_density(im,Tm,MAT);
% pass density to nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom);
rho  = n2interp(1).data;

disp('- INITIALIZING DIKE ARRAYS -');

dm=zeros(size(xm));

% Update dike condition
[BC,ep,epsIIm,epNH,sxxm,Tm,T,im,p]=SiStER_update_dike(BC,X,Y,x,y,xc,yc,T,topo_x,topo_y,xm,ym,im,xsize,ysize,icn,jcn,qd,Tm,ep,epsIIm,epNH,sxxm,0,p,MAT,BCtherm,0,dt_m);


if PARAMS.BalanceStickyLayer~=0 && BC.DIKE.on==1
	PARAMS.BalanceStickyLayer=0; % compensated inflow is incompatible with diking problem
	disp('--------------------------------------------');	
	disp('Balancing sticky air flag incompatible with ');
	disp('internal diking condition - do not worry, ');
	disp('I fixed this for you.');
	disp('--------------------------------------------');
end


