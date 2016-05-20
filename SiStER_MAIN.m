% SiStER_MAIN.m
%
% Simple Stokes solver with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
% jaolive <at> ldeo.columbia.edu
% March 2011 - May 2016

clear all
close all

% INITIALIZATION

% Input File: loads parameter values, model geometry, boundary conditions
InpFil = input('Input file ? ','s');
run(InpFil)

% construct grid and marker arrays
SiStER_Initialize

% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0;
numpicard=zeros(1,Nt); % to keep track of the number of Picard iterations over time



for t=1:Nt % time loop
    

    % Here we prepare nodal arrays to feed the solver and get a (vx,vy,P)
    % solution
    
    % GET MARKER DENSITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % on markers
    [rhom]=SiStER_get_density(im,Tm,MAT);
    % pass density to shear nodes
    [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom);
    rho(:,:)  = n2interp(1).data;

    % GET MARKER ELASTIC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % on markers
    [Gm]=SiStER_get_elastic_moduli(im,MAT);
    % pass shear modulus to normal and shear nodes
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,1./Gm);
    Gn(:,:)=1./(n2interp(1).data);
    [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,1./Gm);
    Gs(:,:) = 1./(n2interp(1).data);   

    % PASS OLD STRESSES TO NODES
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
    sxxOLD(:,:) = n2interp(1).data;
    [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,sxym);
    sxyOLD(:,:) = n2interp(1).data;  

    
    % special case for the very first time iteration 
    if t==1
        % at the very first time step we do not have a pressure field to feed
        % into a yield criterion, so a lithostatic pressure field is
        % approximated from the density structure.
        dY=repmat([0 dy]',1,Nx);
        Plith=cumsum(rho.*dY)*PARAMS.gy;
        [pm]=SiStER_interp_normal_nodes_to_markers(Plith,xc,yc,xm,ym,icn,jcn);
        % at the very first time step we do not have a velocity field
        % either, so the following fields need to be initialized:
        vold=ones(Ny,Nx);
        EPSIIold=ones(Ny-1,Nx-1);
    end
    
    % PICARD ITERATIONS FOR NON-LINEAR FLOW LAWS
    % these make sure that strain rate and viscosities are consistent
    SiStER_run_Picard_iterations;

    % at this point we should have a (vx,vy,p) solution consistent with 
    % the viscosity and stress fields
       
    % UPDATE ELASTIC AND PLASTIC STRESSES (to allow build-up of elastic stresses)
    if (t > 1 && (PARAMS.YNElast==1  || PARAMS.YNPlas==1))
        SiStER_stress_update;
    end
 
    % SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);

    % THERMAL SOLVE  
    if PARAMS.Tsolve==1
        SiStER_thermal_update;
    end
  
    % OUTPUT VARIABLES OF INTEREST (prior to advection)
    if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt % SAVING SELECTED OUTPUT
        disp('SAVING SELECTED VARIABLES') 
        filename=num2str(t);
        save(filename,'X','Y','vx','vy','p','time','xm','ym','etam','rhom','BC','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','syield','pm','icn','jcn','qd','topo_x','topo_y')
    end
    
     
    % MARKER UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SiStER_update_markers;
    % advect markers
    % remove markers if necessary
    % add markers if necessary
    SiStER_update_topography_markers
    % here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get current nodal strain rate on updated markers
    [exxm]=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn);
    [exym]=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn);
    epsIIm=sqrt(exxm.^2+exym.^2);
    % get current pressure field on updated markers
    [pm]=SiStER_interp_normal_nodes_to_markers(p,xc,yc,xm,ym,icn,jcn);

    % update time
    time=time+dt_m;

    disp('---------------')
    disp(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
    disp('--------------------------------')
    disp('--------------------------------')
    

end

disp('FIN')

    