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



% all this is in the update function now - TMorrow 26 Sep 2019

% set up initial dike location, full length of elastic lithosphere
%BC.DIKE.M(X<=BC.DIKE.xR & X>=BC.DIKE.xL)=1;

% M array 0 above Dshallow
%BC.DIKE.M(Y<=(mean(topo_y)+BC.DIKE.Dshallow))=0;

% M array hotter than dike bounding T
%BC.DIKE.M(T>=BC.DIKE.Tmax)=0;

% Set constant M value in dike for now, rated to half the extension rate
%BC.DIKE.M=BC.DIKE.M*BC.DIKE.mval*(abs(BC.right(3))+abs(BC.left(3)))/2;

% Initialize thermal profile
%A=0.006819*(BC.DIKE.Lmin/1e3)^(2.05)+0.2356;
%Tm((abs(xm-xsize/2)<BC.DIKE.Nwide))=1300.*erf(abs(ym(abs(xm-xsize/2)<BC.DIKE.Nwide)-mean(topo_y))./2./sqrt((2*MAT(3).k/MAT(3).rho0/MAT(3).cp).*(A*2.9e6*24*3600*365.25+abs((xm(abs(xm-xsize/2)<BC.DIKE.Nwide))-xsize/2)/BC.right(3))));
%Tm((abs(xm-xsize/2)>=BC.DIKE.Nwide))=1300.*erf(abs(ym((abs(xm-xsize/2)>=BC.DIKE.Nwide))-mean(topo_y))./2./sqrt((2*MAT(3).k/MAT(3).rho0/MAT(3).cp).*(A*2.9e6*24*3600*365.25+abs(BC.DIKE.Nwide)/BC.right(3))));
%Tm(im==1)=PARAMS.T0;

%[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
%T(:,:)  = n2interp(1).data;

%disp('T initialized')
