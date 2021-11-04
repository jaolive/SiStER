function[BC,dm,ep,epsIIm,epNH,sxxm,Tm,T,im,mcreate,p,sxxOLD]=SiStER_update_dike(BC,X,Y,x,y,xc,yc,T,topo_x,topo_y,xm,ym,im,xsize,ysize,icn,jcn,qd,dm,Tm,ep,epsIIm,epNH,sxxm,sxxOLD,t,p,MAT,BCtherm,time,mcreate,dt_m,rho,pp)

% Updates dike parameters (or initializes) each timestep
% TMorrow 26 Sep 2019

% dikes only work in the domain (away from edges) and within a fixed grid size region. If a dike is adjacent to the varaible grid size change, there will be problems

% ------- find topo interface, calculate compensating inflow for outflowing air and rock layers

bL=interp1(topo_x,topo_y,0);
bR=interp1(topo_x,topo_y,xsize);
utop=BC.right(3)*(bL+bR)/xsize;
ubot=BC.right(3)*(2*ysize-bL-bR)/xsize;

% ------- force thermal condition

% shift profile down to match lithosphere base
%A=1.175e-8*(BC.DIKE.Lmin/1e3)^(5.511)+0.9962;
%A=7.616e-13*(BC.DIKE.Lmin/1e3)^(7.696)+0.9962;
%A=0.006819*(BC.DIKE.Lmin/1e3)^(2.05)+0.2356;

% ------- old cooling curve
%Tm=1250.*erf(abs(ym)-mean(topo_y))./2./sqrt((2*MAT(3).k/MAT(3).rho0/MAT(3).cp).*(1.0e6*24*3600*365.25+abs((xm))-xsize/2)/BC.right(3));
%Tm((abs(xm-xsize/2)<BC.DIKE.Nwide))=1300.*erf(abs(ym(abs(xm-xsize/2)<BC.DIKE.Nwide)-mean(topo_y))./2./sqrt((2*MAT(3).k/MAT(3).rho0/MAT(3).cp).*(A*2.9e6*24*3600*365.25+abs((xm(abs(xm-xsize/2)<BC.DIKE.Nwide))-xsize/2)/BC.right(3))));
%Tm((abs(xm-xsize/2)>=BC.DIKE.Nwide))=1300.*erf(abs(ym((abs(xm-xsize/2)>=BC.DIKE.Nwide))-mean(topo_y))./2./sqrt((2*MAT(3).k/MAT(3).rho0/MAT(3).cp).*(A*2.9e6*24*3600*365.25+abs(BC.DIKE.Nwide)/BC.right(3))));

% ------- notch parameters with thermal control
if BC.DIKE.Notch==1

	% T(x,z)=z*700*((1/Lmax-1/Lmin)*x/x_notch+1/Lmin)
	Tm(abs(xm-xsize/2) < BC.DIKE.Nwide)=(ym(abs(xm-xsize/2) < BC.DIKE.Nwide)-mean(topo_y)).*700.*((1/BC.DIKE.Lmax-1/BC.DIKE.Lmin).*abs(xm(abs(xm-xsize/2) < BC.DIKE.Nwide)-xsize/2)./BC.DIKE.Nwide+1/BC.DIKE.Lmin);
	
	Tm(abs(xm-xsize/2) >= BC.DIKE.Nwide)=(ym(abs(xm-xsize/2) >= BC.DIKE.Nwide)-mean(topo_y)).*700.*((1/BC.DIKE.Lmax-1/BC.DIKE.Lmin)+1/BC.DIKE.Lmin);
	

end

% bound temperature
Tm(Tm>BCtherm.bot(2))=BCtherm.bot(2);
Tm(Tm<0)=0;

% get temperature on nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
T(:,:)  = n2interp(1).data;

% check dike timing
if time>=BC.DIKE.tstart*1e6*365.25*24*3600 && time>BC.DIKE.tdynamic*1e6*365.25*24*3600
    
    % ------- update dike arrays
        %disp('-------------^-----------')
 		disp('-  UPDATE - /|\ - DIKE  -')
     	%disp('----------- |~| ---------')

    % update BC.DIKE.M array
    %BC.DIKE.M(X<=BC.DIKE.xR & X>=BC.DIKE.xL)=1;
    
    %BC.DIKE.Tmin=0*BC.DIKE.M+1e20.*(BC.DIKE.M==0);

    %DWColumn=BC.DIKE.DX(find(BC.DIKE.M>0));
    %dikewidth=DWColumn(1);

    % zero M array greater than dike height
    %BC.DIKE.M(Y<=mean(topo_y(topo_x<BC.DIKE.xR & topo_x>BC.DIKE.xL))+BC.DIKE.Dshallow)=0;
    %BC.DIKE.M(
    %BC.DIKE.M(Y<=(mean(topo_y(topo_x<(BC.DIKE.xL+dikewidth) & topo_x>BC.DIKE.xL))+BC.DIKE.Dshallow))=0;

        % initialize
        BC.DIKE.H=[];
        
        
    for ii=2:length(X(1,:))-1
        BC.DIKE.H(1,ii)=interp1(topo_x,topo_y,X(1,ii));%mean(topo_y(topo_x>=X(1,ii) & topo_x<=(X(1,ii+1))));
    end

        BC.DIKE.H(length(BC.DIKE.H)+1)=BC.DIKE.H(end);
        BC.DIKE.H(1)=BC.DIKE.H(2);
        BC.DIKE.H=Y-repmat(BC.DIKE.H,length(Y(:,1)),1);
        %BC.DIKE.H(BC.DIKE.H<0)=0;

    % zero M array hotter than dike bounding T
    %BC.DIKE.M(T>=BC.DIKE.Tmax)=0;
    %HTcut=min(BC.DIKE.H(T>=BC.DIKE.Tmax));
    %BC.DIKE.M(BC.DIKE.M>0 & BC.DIKE.H>(HTcut))=0;
    
    %BC.DIKE.M((max(ind)+1):(max(ind)+10))=1;
    
    % ------- Testing a two element dike
    
    %BC.DIKE.M(X<BC.DIKE.xR & X>BC.DIKE.xL & BC.DIKE.H>0 & T<800)=0.5/2;
    
    %BC.DIKE.MV=BC.DIKE.M*(abs(BC.right(3)))*2; % get M*V
    
    
    % ------- Iterative guess and check
%     
%     BC.DIKE.dsxxSURR=conv2(BC.DIKE.dsxxPic,BC.DIKE.KER,'same');
%     
     BC.DIKE.SXX=BC.DIKE.sxxOLDPic;
     BC.DIKE.CONVSXX=conv2(BC.DIKE.sxxOLDPic,BC.DIKE.KER,'same');
     BC.DIKE.SXX(:,161)=BC.DIKE.CONVSXX(:,161);
 
     BC.DIKE.P=MAT(2).rho0*9.8.*BC.DIKE.H; %=2800;p; % lithostatic pressure
     
     BC.DIKE.Pm=BC.DIKE.rho*9.8.*BC.DIKE.H+BC.DIKE.P0; % magma pressure
     
%     BC.DIKE.dsxxMOLD=BC.DIKE.dsxxM;
     
     BC.DIKE.dsxxM=BC.DIKE.SXX-BC.DIKE.P+BC.DIKE.Pm-BC.DIKE.T; % tensile stress excess
%     %BC.DIKE.dsxxM(BC.DIKE.dsxxM<-1e12)=0;
    
     if pp==1
         %BC.DIKE.M(BC.DIKE.dsxxM>0)=0.5;
         BC.DIKE.sxxLIM=1.2e9; %max(max(BC.DIKE.dsxxM)); % max tensile stress excess at the start of the solve controls our initial guess
         %BC.DIKE.DelM=BC.DIKE.dsxxM./max(max(BC.DIKE.dsxxM)); % proportion decides initial M guess
         BC.DIKE.DelM(BC.DIKE.H>0 & T<1000)=(rand(size(BC.DIKE.DelM(BC.DIKE.H>0 & T<1000)))+0.5)*0.25;%0.25;
         BC.DIKE.DelM(X>BC.DIKE.xR)=0;
         BC.DIKE.DelM(X<BC.DIKE.xL)=0;
         BC.DIKE.DelM(BC.DIKE.H<0)=0;
%     %elseif mod(pp,2)>0
%     %    BC.DIKE.M(BC.DIKE.dsxxM>0)=BC.DIKE.Mold(BC.DIKE.dsxxM>0)+0.1;
%     %    BC.DIKE.M(BC.DIKE.dsxxM<0)=BC.DIKE.Mold(BC.DIKE.dsxxM<0)-0.1;
%     %    BC.DIKE.M(BC.DIKE.M>1)=1;
%     %    BC.DIKE.M(BC.DIKE.M<0)=0;
%         %BC.DIKE.DelM=-(BC.DIKE.dsxxM-BC.DIKE.dsxxMOLD)./2./BC.DIKE.sxxLIM; % determine del M
%         %BC.DIKE.DelM(BC.DIKE.T>1e19)=0; % reinforce the BC.DIKE.T location limit
%         %BC.DIKE.DelM(BC.DIKE.H>0 & T<800 & BC.DIKE.T<1e9)=0.1;
     else
         %BC.DIKE.dsxxSURR=conv2(BC.DIKE.dsxxPic,BC.DIKE.KER,'same');
         BC.DIKE.dsxxRAT=BC.DIKE.dsxxM./BC.DIKE.sxxLIM;
         BC.DIKE.DelM(:)=0;
         BC.DIKE.DelM=BC.DIKE.dsxxRAT;
         BC.DIKE.DelM(X>BC.DIKE.xR)=0;
         BC.DIKE.DelM(X<BC.DIKE.xL)=0;
         BC.DIKE.DelM(BC.DIKE.H<0)=0;
% 
%         BC.DIKE.M(BC.DIKE.Mold>0)=BC.DIKE.Mold(BC.DIKE.Mold>0)-BC.DIKE.Mold(BC.DIKE.Mold>0).*BC.DIKE.dsxxSURR(BC.DIKE.Mold>0) % ./min(BC.DIKE.dsxxSURR(BC.DIKE.Mold>0));
%         
%         BC.DIKE.DelM=BC.DIKE.M-BC.DIKE.Mold;
%         
%         %BC.DIKE.DelM=-(BC.DIKE.dsxxM-BC.DIKE.dsxxMOLD)./2./BC.DIKE.sxxLIM; % determine del M
%         %BC.DIKE.DelM(X>BC.DIKE.xR)=0;
%         %BC.DIKE.DelM(X<BC.DIKE.xL)=0;
%         %BC.DIKE.DelM(BC.DIKE.H<0)=0;
%         %BC.DIKE.DelM(BC.DIKE.H>0 & T<600 & T>300 & BC.DIKE.T<1e9)=(-0.5+rand(size(BC.DIKE.DelM(BC.DIKE.H>0 & T<600 & T>300 & BC.DIKE.T<1e9))))*0.2;
     end
%    
     BC.DIKE.MOld=BC.DIKE.M;
     BC.DIKE.M=BC.DIKE.M+BC.DIKE.DelM;
% 
%     %BC.DIKE.M=BC.DIKE.dsxxM./max(max(BC.DIKE.dsxxM)); % scale M to max tensile stress excess
%     
%     %BC.DIKE.M=BC.DIKE.M.*0.1; % max dike injection is 0.1
%     
     BC.DIKE.M(BC.DIKE.M<0)=0; % zero non-injecting nodes
% 
%     %BC.DIKE.M=BC.DIKE.Mold+0.2*BC.DIKE.DelM./max(max(BC.DIKE.dsxxM)); % add del M
%     %                       ^ damping factor
% 
     BC.DIKE.M(BC.DIKE.M>1)=1; % cap M at 1
     
     BC.DIKE.M=movmean(BC.DIKE.M,BC.DIKE.SmoothCells); % smoothing we like, TMorrow 02 Sept. 2021
% 
     BC.DIKE.M(X>BC.DIKE.xR)=0;
     BC.DIKE.M(X<BC.DIKE.xL)=0;
     BC.DIKE.M(BC.DIKE.H<0)=0; % zero air injection, TMorrow 07 Oct 2021
%      
%     %BC.DIKE.M(BC.DIKE.M<0)=0; % zero non-injecting nodes
%     
     BC.DIKE.MV=BC.DIKE.M*(abs(BC.right(3))+abs(BC.left(3))); % get M*V
%     
%     BC.DIKE.OPEN(BC.DIKE.MV>0)=1; % dike opening flag (helpful information but otherwise unused?)
%     
%     %BC.DIKE.DelM=BC.DIKE.M-BC.DIKE.Mold; % visualize change in M with iteration
%     
%     BC.DIKE.Mold=BC.DIKE.M; % hold on to old M
% 
%     
%     
%     %sxxOLD=sxxOLD-0.2*BC.DIKE.DelM; % correct sxx
%     %              ^ damping factor again

    % ---------------------
    

    % ------- dike pressure
    %BC.DIKE.Pm=BC.DIKE.rho*9.8.*BC.DIKE.H-BC.DIKE.P0*1e6;

    %BC.DIKE.Pm(BC.DIKE.rho*9.8.*BC.DIKE.H<0)=0;

    %airdike=find(BC.DIKE.H<0 & BC.DIKE.M>0);

    % taper M at top of dike
    %ind=find(BC.DIKE.M>0);
    %dtop=airdike(1);

    %BC.DIKE.M(dtop-5)=0.05;
    %BC.DIKE.M(dtop-4)=0.25;
    %BC.DIKE.M(dtop-3)=0.5;
    %BC.DIKE.M(dtop-2)=0.75;
    %BC.DIKE.M(dtop-1)=1.0;

	%BC.DIKE.SXX=sxxOLD;

    %ind=find(BC.DIKE.M);

    %if t>2 && BC.DIKE.P0>0
        

%         if t>3
%             BC.DIKE.SXX=(BC.DIKE.SXX+sxxOLD)./2;
%         else
%             BC.DIKE.SXX=sxxOLD;            
%         end
%            
        %BC.DIKE.SXX=sxxOLD;
        %BC.DIKE.SXX=conv2(sxxOLD,BC.DIKE.KER,'same');
        %BC.DIKE.P=p;
        %BC.DIKE.P=cumsum(rho*9.8.*BC.DIKE.DY);
		%BC.DIKE.RhoMain=rho;
        %BC.DIKE.Pm=BC.DIKE.rho*9.8.*BC.DIKE.H+BC.DIKE.P0;
        
        % ------- Diffuse yielding setup
        
        %BC.DIKE.K=2*MAT(4).G*(1+0.25)/(3*(1-2*0.25))/100;
        %BC.DIKE.dsxxMold=BC.DIKE.dsxxM;
        %BC.DIKE.dsxxM=BC.DIKE.SXX-BC.DIKE.P+BC.DIKE.Pm-BC.DIKE.T; % dsxxM excess
        %BC.DIKE.dsxxM(1,:)=0; % model boundaries get no correction, always
        %BC.DIKE.dsxxM(:,end)=0; % model boundaries get no correction, always
        
        % ---- smoothing
        %KER=[0.5 1 1 1 1 1 0.5]'./6;
        %BC.DIKE.dsxxM=conv2(BC.DIKE.dsxxM,BC.DIKE.KER,'same');

	%BC.DIKE.dsxxM(T>BC.DIKE.Tmax)=0;
        
        
        
        
        %BC.DIKE.DELdsxxM=0.2*(BC.DIKE.dsxxM-BC.DIKE.OLDdsxxM);
        %BC.DIKE.DELdsxxM(BC.DIKE.DELdsxxM<0)=0;
        
        %BC.DIKE.dsxxM=BC.DIKE.dsxxM./5;   
        
        % ---- vertical smoothing test 29 June 2021
        %BC.DIKE.dsxxM=smoothdata(BC.DIKE.dsxxM,'Gaussian',20);
        
        
        % ------- set M from magmatic elasticity and stress excess
        %sxxOLD=sxxOLD-BC.DIKE.dsxxM;
        %sxxOLD(sxxOLD<0 & BC.DIKE.dsxxM>0)=0;
        
        %BC.DIKE.EpM=(BC.DIKE.dsxxM)/BC.DIKE.K/dt_m;
               
        %BC.DIKE.M=BC.DIKE.EpM.*BC.DIKE.DX/(2*BC.right(3));
        
        

        %sxxOLD=sxxOLD-BC.DIKE.dsxxM;
        
        
        %sxxOLD(BC.DIKE.dsxxM>0)=0;
        %sxxOLD(sxxOLD<0 & BC.DIKE.dsxxM>0)=0;
        
        
        
        
        %BC.DIKE.EpM=(BC.DIKE.DELdsxxM)/BC.DIKE.K/dt_m;
        %BC.DIKE.M=BC.DIKE.EpM.*BC.DIKE.DX/(2*BC.right(3))+BC.DIKE.M;
        
        %sxxOLD=sxxOLD+BC.DIKE.DELdsxxM;
        
        %BC.DIKE.OLDdsxxM=conv2(sxxOLD,BC.DIKE.KER,'same');
               
        %BC.DIKE.dM=BC.DIKE.EpM.*BC.DIKE.DX/(2*BC.right(3))-BC.DIKE.M;
        
        %BC.DIKE.M=BC.DIKE.M+0.2*BC.DIKE.dM;
        
        %BC.DIKE.M(BC.DIKE.M>1)=1.0;
        
        %BC.DIKE.M(:)=0;
        
    % ----- imposing M to watch stress perturbation 19 July 2021
        % TMorrow
        %BC.DIKE.M(50:55,161)=1.0;
        
    %if BC.DikeKillFlag==0
    %    BC.DIKE.M(X<BC.DIKE.xR & X>BC.DIKE.xL & BC.DIKE.H>0 & T<800)=BC.DIKE.M(X<BC.DIKE.xR & X>BC.DIKE.xL & BC.DIKE.H>0 & T<800)+0.025;
    %    
    %    if mean(BC.DIKE.M(BC.DIKE.M>0))>1
    %        BC.DikeKillFlag=1;
    %        BC.DIKE.M(:)=0;
    %    end
    %end
    
    %BC.DIKE.MV=BC.DIKE.M*(abs(BC.right(3)))*2; %BC.DIKE.MV is the dike RHS term
        
        % temporary fixed M to test dilatational strain
	%BC.DIKE.M(BC.DIKE.H>0 & T<800 & BC.DIKE.T<1e9)=0.8; 
	%BC.DIKE.M(T<BC.DIKE.Tmax & BC.DIKE.T<1e9 & BC.DIKE.H>=0)=BC.DIKE.mval;     % dike through domain
    
    % scale dike for multiple columns
    %BC.DIKE.SCAL=repmat(sum(BC.DIKE.M,2),1,numel(BC.DIKE.SCAL(1,:)));
    %BC.DIKE.M=BC.DIKE.M./BC.DIKE.SCAL;
       
    % ------- update dike markers
    %dm(:)=0;
	%im(dm==1 & im>1)=2;
    
    % ------- heal advected dike material from last time step
%	ep(dm>0)=0;
    
    % ------- redefine dike markers
	dm(:)=0;

    % find markers in/around dike area (allow for interpolation bleed to counter strain interpolation bleed?)
    [dm]=SiStER_interp_normal_nodes_to_markers(BC.DIKE.M,xc,yc,xm,ym,icn,jcn);

    dm(dm>0)=1;
    im(dm==1 & im>1)=4;
    %im(dm>0 & 
    mcreate(dm==1 & mcreate==0)=time;

    %airind=interp1(topo_x,topo_y,(BC.DIKE.xR+BC.DIKE.xL)/2);
    %im(ym<airind & dm>0)=1;
    
    % ---------------------
    % M*V value in dike and opening flag

    %BC.DIKE.M=BC.DIKE.M*BC.DIKE.mval;
    %BC.DIKE.MV=BC.DIKE.M*(abs(BC.right(3))+abs(BC.left(3)));
    %BC.DIKE.MEXX=BC.DIKE.MV./[repmat(dx,length(BC.DIKE.M(:,1)),1) zeros(length(BC.DIKE.M(:,1)),1)];
    %BC.DIKE.OPEN(BC.DIKE.MV>0)=1;

    % ---------------------
	% compensate dike source term in inflow from bottom

	airdike=find(BC.DIKE.H<0 & BC.DIKE.MV>0);

	crustdike=find(BC.DIKE.H>=0 & BC.DIKE.MV>0);

	%  second term is dike height times M value times velocity, in locations where dike exists
	ubot=ubot-sum(BC.DIKE.DY(crustdike).*BC.DIKE.MV(crustdike))/xsize;
	utop=utop-sum(BC.DIKE.DY(airdike).*BC.DIKE.MV(airdike))/xsize;

    % ---------------------
    
    % ------- Adjust temperature bounds
    %Tm(Tm>BCtherm.bot(2))=BCtherm.bot(2);
    %Tm(im==1)=0;
    
    % --------------------- 
    
    % ------- get temperature on nodes
    %[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
    %T(:,:)  = n2interp(1).data;

    % ---------------------

    % Add dike temperature
    [dmv]=SiStER_interp_shear_nodes_to_markers(BC.DIKE.MV,x,y,xm,ym,icn,jcn);
    [dmx]=SiStER_interp_shear_nodes_to_markers(BC.DIKE.DX,x,y,xm,ym,icn,jcn);
    %Tsol < Tdike < Tliq
    %Tm(dmv>0 & Tm>=1000 & Tm<1300) = Tm(dmv>0 & Tm>=1000 & Tm<1300)+dmv(dmv>0 & Tm>=1000 & Tm<1300)./dmx(dmv>0 & Tm>=1000 & Tm<1300).*(1300-Tm(dmv>0 & Tm>=1000 & Tm<1300)).*dt_m;
    %Tdike < Tsol
    %Tm(dmv>0 & Tm<1000) = Tm(dmv>0 & Tm<1000)+dmv(dmv>0 & Tm<1000)./dmx(dmv>0 & Tm<1000).*(1300-Tm(dmv>0 & Tm<1000) + 5e5/1e3).*dt_m;

    % ---------------------

    % ------- heal dike
	%ep(dm>0)=0;

    % ---------------------
    
    % ------- debug plots
    %plot(BC.DIKE.M(:,161))
    % drawnow
    % hold on

% elseif time>=BC.DIKE.tstart*1e6*365.25*24*3600
%     
%     BC.DIKE.DelM(BC.DIKE.H>0 & T<800 & BC.DIKE.T<1e9)=0.1;
%     
%     BC.DIKE.M=BC.DIKE.DelM;
% 
%     %BC.DIKE.M=BC.DIKE.dsxxM./max(max(BC.DIKE.dsxxM)); % scale M to max tensile stress excess
%     
%     %BC.DIKE.M=BC.DIKE.M.*0.1; % max dike injection is 0.1
%     
%     BC.DIKE.M(BC.DIKE.M<0)=0; % zero non-injecting nodes
% 
%     %BC.DIKE.M=BC.DIKE.Mold+0.2*BC.DIKE.DelM./max(max(BC.DIKE.dsxxM)); % add del M
%     %                       ^ damping factor
% 
%     BC.DIKE.M(BC.DIKE.M>1)=1; % cap M at 1
% 
%     BC.DIKE.M(X>BC.DIKE.xR)=0;
%     BC.DIKE.M(X<BC.DIKE.xL)=0;
%      
%     %BC.DIKE.M(BC.DIKE.M<0)=0; % zero non-injecting nodes
%     
%     BC.DIKE.MV=BC.DIKE.M*(abs(BC.right(3))+abs(BC.left(3))); % get M*V
%     
%     BC.DIKE.OPEN(BC.DIKE.MV>0)=1; % dike opening flag (helpful information but otherwise unused?)
%     
%     %BC.DIKE.DelM=BC.DIKE.M-BC.DIKE.Mold; % visualize change in M with iteration
%     
%     BC.DIKE.Mold=BC.DIKE.M; % hold on to old M
% 
%     
%     
%     %sxxOLD=sxxOLD-0.2*BC.DIKE.DelM; % correct sxx
%     %              ^ damping factor again
% 
%     % ---------------------
%     
% 
%     % ------- dike pressure
%     %BC.DIKE.Pm=BC.DIKE.rho*9.8.*BC.DIKE.H-BC.DIKE.P0*1e6;
% 
%     %BC.DIKE.Pm(BC.DIKE.rho*9.8.*BC.DIKE.H<0)=0;
% 
%     %airdike=find(BC.DIKE.H<0 & BC.DIKE.M>0);
% 
%     % taper M at top of dike
%     %ind=find(BC.DIKE.M>0);
%     %dtop=airdike(1);
% 
%     %BC.DIKE.M(dtop-5)=0.05;
%     %BC.DIKE.M(dtop-4)=0.25;
%     %BC.DIKE.M(dtop-3)=0.5;
%     %BC.DIKE.M(dtop-2)=0.75;
%     %BC.DIKE.M(dtop-1)=1.0;
% 
% 	%BC.DIKE.SXX=sxxOLD;
% 
%     %ind=find(BC.DIKE.M);
% 
%     %if t>2 && BC.DIKE.P0>0
%         
% 
% %         if t>3
% %             BC.DIKE.SXX=(BC.DIKE.SXX+sxxOLD)./2;
% %         else
% %             BC.DIKE.SXX=sxxOLD;            
% %         end
% %            
%         %BC.DIKE.SXX=sxxOLD;
%         %BC.DIKE.SXX=conv2(sxxOLD,BC.DIKE.KER,'same');
%         %BC.DIKE.P=p;
%         %BC.DIKE.P=cumsum(rho*9.8.*BC.DIKE.DY);
% 		%BC.DIKE.RhoMain=rho;
%         %BC.DIKE.Pm=BC.DIKE.rho*9.8.*BC.DIKE.H+BC.DIKE.P0;
%         
%         % ------- Diffuse yielding setup
%         
%         %BC.DIKE.K=2*MAT(4).G*(1+0.25)/(3*(1-2*0.25))/100;
%         %BC.DIKE.dsxxMold=BC.DIKE.dsxxM;
%         %BC.DIKE.dsxxM=BC.DIKE.SXX-BC.DIKE.P+BC.DIKE.Pm-BC.DIKE.T; % dsxxM excess
%         %BC.DIKE.dsxxM(1,:)=0; % model boundaries get no correction, always
%         %BC.DIKE.dsxxM(:,end)=0; % model boundaries get no correction, always
%         
%         % ---- smoothing
%         %KER=[0.5 1 1 1 1 1 0.5]'./6;
%         %BC.DIKE.dsxxM=conv2(BC.DIKE.dsxxM,BC.DIKE.KER,'same');
% 
% 	%BC.DIKE.dsxxM(T>BC.DIKE.Tmax)=0;
%         
%         
%         
%         
%         %BC.DIKE.DELdsxxM=0.2*(BC.DIKE.dsxxM-BC.DIKE.OLDdsxxM);
%         %BC.DIKE.DELdsxxM(BC.DIKE.DELdsxxM<0)=0;
%         
%         %BC.DIKE.dsxxM=BC.DIKE.dsxxM./5;   
%         
%         % ---- vertical smoothing test 29 June 2021
%         %BC.DIKE.dsxxM=smoothdata(BC.DIKE.dsxxM,'Gaussian',20);
%         
%         
%         % ------- set M from magmatic elasticity and stress excess
%         %sxxOLD=sxxOLD-BC.DIKE.dsxxM;
%         %sxxOLD(sxxOLD<0 & BC.DIKE.dsxxM>0)=0;
%         
%         %BC.DIKE.EpM=(BC.DIKE.dsxxM)/BC.DIKE.K/dt_m;
%                
%         %BC.DIKE.M=BC.DIKE.EpM.*BC.DIKE.DX/(2*BC.right(3));
%         
%         
% 
%         %sxxOLD=sxxOLD-BC.DIKE.dsxxM;
%         
%         
%         %sxxOLD(BC.DIKE.dsxxM>0)=0;
%         %sxxOLD(sxxOLD<0 & BC.DIKE.dsxxM>0)=0;
%         
%         
%         
%         
%         %BC.DIKE.EpM=(BC.DIKE.DELdsxxM)/BC.DIKE.K/dt_m;
%         %BC.DIKE.M=BC.DIKE.EpM.*BC.DIKE.DX/(2*BC.right(3))+BC.DIKE.M;
%         
%         %sxxOLD=sxxOLD+BC.DIKE.DELdsxxM;
%         
%         %BC.DIKE.OLDdsxxM=conv2(sxxOLD,BC.DIKE.KER,'same');
%                
%         %BC.DIKE.dM=BC.DIKE.EpM.*BC.DIKE.DX/(2*BC.right(3))-BC.DIKE.M;
%         
%         %BC.DIKE.M=BC.DIKE.M+0.2*BC.DIKE.dM;
%         
%         %BC.DIKE.M(BC.DIKE.M>1)=1.0;
%         
%         %BC.DIKE.M(:)=0;
%         
%         % ----- imposing M to watch stress perturbation 19 July 2021
%         % TMorrow
%         %BC.DIKE.M(50:55,161)=1.0;
%         
%         % temporary fixed M to test dilatational strain
% 	%BC.DIKE.M(BC.DIKE.H>0 & T<800 & BC.DIKE.T<1e9)=0.8; 
% 	%BC.DIKE.M(T<BC.DIKE.Tmax & BC.DIKE.T<1e9 & BC.DIKE.H>=0)=BC.DIKE.mval;     % dike through domain
%     
%     % scale dike for multiple columns
%     %BC.DIKE.SCAL=repmat(sum(BC.DIKE.M,2),1,numel(BC.DIKE.SCAL(1,:)));
%     %BC.DIKE.M=BC.DIKE.M./BC.DIKE.SCAL;
%        
%     % ------- update dike markers
%     %dm(:)=0;
% 	im(dm==1 & im>1)=2;
%     
%     % ------- heal advected dike material from last time step
% 	ep(dm>0)=0;
%     
%     % ------- redefine dike markers
% 	dm(:)=0;
% 
%     % find markers in/around dike area (allow for interpolation bleed to counter strain interpolation bleed?)
%     [dm]=SiStER_interp_normal_nodes_to_markers(BC.DIKE.M,xc,yc,xm,ym,icn,jcn);
% 
%     dm(dm>0)=1;
%     im(dm==1 & im>1)=4;
% 	mcreate(dm==1 & mcreate==0)=time;
% 
%     %airind=interp1(topo_x,topo_y,(BC.DIKE.xR+BC.DIKE.xL)/2);
%     %im(ym<airind & dm>0)=1;
%     
%     % ---------------------
%     % M*V value in dike and opening flag
% 
%     %BC.DIKE.M=BC.DIKE.M*BC.DIKE.mval;
%     %BC.DIKE.MV=BC.DIKE.M*(abs(BC.right(3))+abs(BC.left(3)));
%     %BC.DIKE.MEXX=BC.DIKE.MV./[repmat(dx,length(BC.DIKE.M(:,1)),1) zeros(length(BC.DIKE.M(:,1)),1)];
%     %BC.DIKE.OPEN(BC.DIKE.MV>0)=1;
% 
%     % ---------------------
% 	% compensate dike source term in inflow from bottom
% 
% 	airdike=find(BC.DIKE.H<0 & BC.DIKE.MV>0);
% 
% 	crustdike=find(BC.DIKE.H>=0 & BC.DIKE.MV>0);
% 
% 	%  second term is dike height times M value times velocity, in locations where dike exists
% 	ubot=ubot-sum(BC.DIKE.DY(crustdike).*BC.DIKE.MV(crustdike))/xsize;
% 	utop=utop-sum(BC.DIKE.DY(airdike).*BC.DIKE.MV(airdike))/xsize;
% 
%     % ---------------------
%     
%     % ------- Adjust temperature bounds
%     %Tm(Tm>BCtherm.bot(2))=BCtherm.bot(2);
%     %Tm(im==1)=0;
%     
%     % --------------------- 
%     
%     % ------- get temperature on nodes
%     %[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
%     %T(:,:)  = n2interp(1).data;
% 
%     % ---------------------
% 
%     % Add dike temperature
%     [dmv]=SiStER_interp_shear_nodes_to_markers(BC.DIKE.MV,x,y,xm,ym,icn,jcn);
%     [dmx]=SiStER_interp_shear_nodes_to_markers(BC.DIKE.DX,x,y,xm,ym,icn,jcn);
%     %Tsol < Tdike < Tliq
%     %Tm(dmv>0 & Tm>=1000 & Tm<1300) = Tm(dmv>0 & Tm>=1000 & Tm<1300)+dmv(dmv>0 & Tm>=1000 & Tm<1300)./dmx(dmv>0 & Tm>=1000 & Tm<1300).*(1300-Tm(dmv>0 & Tm>=1000 & Tm<1300)).*dt_m;
%     %Tdike < Tsol
%     %Tm(dmv>0 & Tm<1000) = Tm(dmv>0 & Tm<1000)+dmv(dmv>0 & Tm<1000)./dmx(dmv>0 & Tm<1000).*(1300-Tm(dmv>0 & Tm<1000) + 5e5/1e3).*dt_m;
% 
%     % ---------------------
% 
%     % ------- heal dike
% 	%ep(dm>0)=0;
% 
%     % ---------------------
%     
%     % ------- debug plots
%     %plot(BC.DIKE.M(:,161))
%     % drawnow
%     % hold on
%     
%     
else
    %disp('---------------------------');
    disp('-       - NO DIKE -       -');
    %disp('---------------------------');

        BC.DIKE.OPEN(:)=0;
		BC.DIKE.M(:)=0;
        BC.DIKE.MV(:)=0;
		%BC.DIKE.M(:)=0;
end

% TODO ask BEN about vectorizing this operation
% T Morrow 10 Oct 2019
%dm(icn==(find(X(1,:)==unique(X(find(BC.DIKE.M>0))))) & jcn==(find(Y(:,1)==unique(Y(find(BC.DIKE.M>0)))));
% find markers within dike

%dm(xm<=BC.DIKE.xR & xm>=BC.DIKE.xL & ym>(mean(topo_y(topo_x<BC.DIKE.xR & topo_x>BC.DIKE.xL))+BC.DIKE.Dshallow) & Tm<600)=1;

% ------- force lithosphere mantle

if BC.DIKE.Notch==1
    %im(im==7 & ym<90e3)=6; % fixed - adjust if model domain changes TMorrow 14 Apr 2020
    
    im(im==3 & Tm<=700)=2;
    im(im==2 & Tm>700)=3;
    im(im==4 & Tm>700)=5;
    im(im==5 & Tm<=700)=4;
    
    
    %im(im==3 & ym<=(abs(xm-xsize/2)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)<BC.DIKE.Nwide)=2;% fixed sloped lithosphere
    %im(im==3 & ym<=(abs(BC.DIKE.Nwide)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)>=BC.DIKE.Nwide)=2;% fixed flat lithosphere

    %im(im==2 & ym>(abs(xm-xsize/2)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)<BC.DIKE.Nwide)=3;% fixed sloped lithosphere
    %im(im==2 & ym>(abs(BC.DIKE.Nwide)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)>=BC.DIKE.Nwide)=3;% fixed flat lithosphere

    %im(im==4 & ym>(abs(xm-xsize/2)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)<BC.DIKE.Nwide)=5;% fixed sloped lithosphere
    %im(im==4 & ym>(abs(BC.DIKE.Nwide)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)>=BC.DIKE.Nwide)=5;% fixed flat lithosphere

    %im(im==5 & ym<=(abs(xm-xsize/2)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)<BC.DIKE.Nwide)=4;% fixed sloped lithosphere
    %im(im==5 & ym<=(abs(BC.DIKE.Nwide)/BC.DIKE.Lslope+BC.DIKE.Lmin+mean(topo_y)) & abs(xm-xsize/2)>=BC.DIKE.Nwide)=4;% fixed flat lithosphere
end

% ------- update BC velocities

BC.top(3)=utop;
BC.bot(3)=-ubot;

% ---------------------
