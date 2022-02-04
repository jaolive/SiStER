function[BC,ep,epsIIm,epNH,sxxm,Tm,T,im,p,sxxOLD]=SiStER_update_dike(BC,X,Y,x,y,xc,yc,T,topo_x,topo_y,xm,ym,im,xsize,ysize,icn,jcn,qd,Tm,ep,epsIIm,epNH,sxxm,t,p,MAT,BCtherm,time,dt_m,rho,pp)

% Updates dike parameters (or initializes) each timestep
% TMorrow 26 Sep 2019

% dikes only work in the domain (away from edges) and within a fixed grid size region. If a dike is adjacent to the varaible grid size change, there will be problems

% ------- impose thermal condition if set

%if BC.DIKE.imposethermal>0
%    SiStER_impose_thermal_condition
%end

% ---------------------

% ------- find topo interface, calculate compensating inflow for outflowing air and rock layers

bL=interp1(topo_x,topo_y,0);
bR=interp1(topo_x,topo_y,xsize);
utop=BC.right(3)*(bL+bR)/xsize;
ubot=BC.right(3)*(2*ysize-bL-bR)/xsize;

% ---------------------

% ------- H is depth from surface

BC.DIKE.H=[];

for ii=2:length(X(1,:))-1
    BC.DIKE.H(1,ii)=interp1(topo_x,topo_y,X(1,ii));
end

BC.DIKE.H(length(BC.DIKE.H)+1)=BC.DIKE.H(end);
BC.DIKE.H(1)=BC.DIKE.H(2);
BC.DIKE.H=Y-repmat(BC.DIKE.H,length(Y(:,1)),1);

% ---------------------

% ------- constrain injection location


BC.DIKE.M(:)=BC.DIKE.mval;

BC.DIKE.M(Y<BC.DIKE.top)=0;
BC.DIKE.M(Y>BC.DIKE.bot)=0;
BC.DIKE.M(X>BC.DIKE.xR)=0;
BC.DIKE.M(X<BC.DIKE.xL)=0;

% ---------------------

% ------- M*V is the divergence source term

BC.DIKE.MV=BC.DIKE.M.*(abs(BC.right(3))+abs(BC.left(3))); % get M*V

% ---------------------
    
% ------- track dike material, assign injection phase
dm(:)=0;

% find markers in/around dike area
[dm]=SiStER_interp_normal_nodes_to_markers(BC.DIKE.M,xc,yc,xm,ym,icn,jcn);
dm(dm>0)=1;

% assign injected phase (unless sticky air)
im(dm==1 & im>1)=BC.DIKE.injmat;

%mcreate(dm==1 & mcreate==0)=time;

%airind=interp1(topo_x,topo_y,(BC.DIKE.xR+BC.DIKE.xL)/2);
%im(ym<airind & dm>0)=1;

% ---------------------

% ------- compensate dike source term in inflow from bottom

airdike=find(BC.DIKE.H<0 & BC.DIKE.MV>0);
crustdike=find(BC.DIKE.H>=0 & BC.DIKE.MV>0);

%  second term is dike height times M value times velocity, in locations where dike exists
ubot=ubot-sum(BC.DIKE.DY(crustdike).*BC.DIKE.MV(crustdike))/xsize;
utop=utop-sum(BC.DIKE.DY(airdike).*BC.DIKE.MV(airdike))/xsize;

% ---------------------
    

% ------- update BC velocities

BC.top(3)=utop;
BC.bot(3)=-ubot;

% ---------------------
