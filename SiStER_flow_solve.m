% SiStER FLOW SOLVE
% = one flow solve or one Picard iteration
%
% 1/ calculates viscosity on markers based on current strain rate
% 2/ passes it to nodes, solves Stokes flow
% 3/ extracts strain rate from solution


if PARAMS.BalanceStickyLayer==1
% BALANCE FLUXES %%% JAO July 16, 2015
% RE-ADJUST BCs SO FLUX OF ROCK AND STICKY AIR MATERIAL ARE CONSERVED
% locate height of sticky layer - rock contact on the sides
    bL=interp1(topo_x,topo_y,0);
    bR=interp1(topo_x,topo_y,xsize);
    utop=BC.right(3)*(bL+bR)/xsize;
    ubot=BC.right(3)*(2*ysize-bL-bR)/xsize;
    BC.top(3)=utop;
    BC.bot(3)=-ubot;
end


% GET MARKER VISCOSITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmNonNeg=pm;
pmNonNeg(pmNonNeg<0)=0; % to make sure we don't feed negative pressures in the viscosity
% and yield stress calculation
[etam,syield]=SiStER_get_viscosity(im,Tm,epsIIm,MAT,PARAMS,ep,pmNonNeg);


% at this point the marker densities, viscosities and shear moduli
% are ready to be passed to nodes for the next solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INTERPOLATE VISCOSITY FROM MARKER TO NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,etam);
etan(:,:)=n2interp(1).data;
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,etam);
etas(:,:) = n2interp(1).data;   


% ASSEMBLE SPECIFIC ELASTICITY TERMS (ONLY AFTER FIRST TIME STEP)
if t==1 || PARAMS.YNElast==0
    Zs=ones(size(etas));
    Zn=ones(size(etan));
elseif t>1 && PARAMS.YNElast==1
    Zs=Gs*dt_m./(Gs*dt_m+etas);
    Zn=Gn*dt_m./(Gn*dt_m+etan);
end
% right-hand size (1-Z)*sigmaOLD_ij
srhs_xx=(1-Zn).*sxxOLD;
srhs_xy=(1-Zs).*sxyOLD;


% ASSEMBLE AND SOLVE LINEAR SYSTEM WITH DIRECT SOLVER (BACKSLASH)
[S, Kc, Kb]=SiStER_Stokes_solver(dx,dy,Zs.*etas,Zn.*etan,rho,BC,PARAMS,srhs_xx,srhs_xy); 
[p, vx, vy]=SiStER_reshape_solver_output(S,Kc,Nx,Ny); 
% current solution stored on staggered grid (p, vx, vy)

% extract strain rate from current solution
[EXX,EXY]=SiStER_get_strain_rate(vx,vy,dx,dy,BC); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if PARAMS.YNElast==1
if PARAMS.epsII_from_stress==0 || t==1 % compute strain rate directly from velocity solution
    
    % get current nodal strain rate on (updated) markers
    [exxm]=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn);
    [exym]=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn);
    epsIIm=sqrt(exxm.^2+exym.^2);
    % this strain rate will be used to estimate viscosity for the next Picard iteration
    
elseif PARAMS.epsII_from_stress==1  
    % or compute strain rate directly from the latest stresses
    % to limit the number of interpolations
    % option added by G. Ito, Fall 2015
    % Compute new estimate of stresses (but do not change the saved ones)
    dsxx=(2*etan.*EXX-sxxOLD).*Zn;
    dsxy=(2*etas.*EXY-sxyOLD).*Zs;

    [dsxxm]=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,xm,ym,icn,jcn);
    [dsxym]=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,xm,ym,icn,jcn);
    
    sxxmtemp=sxxm+dsxxm;  
    sxymtemp=sxym+dsxym;   

    % estimate strainrate using updated stresses
    exxm=sxxmtemp./(2*etam)+(dsxxm)./(2*dt_m.*Gm);  
    exym=sxymtemp./(2*etam)+(dsxym)./(2*dt_m.*Gm);
    epsIIm=sqrt(exxm.^2+exym.^2);
    
end
end


