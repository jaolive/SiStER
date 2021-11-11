% SiStER_Input_File


% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
Nt=50; % max number of time iterations
dt_out=1; %10 % output files every "dt_out" iterations


% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize=200e3;
ysize=100e3;
% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y

GRID.dx(1)=1000;
GRID.x(1)=200e3;
GRID.dy(1)=1000;
GRID.y(1)=100e3;

% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=8; % number of markers in the smallest quadrant
Mquad_crit=4; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nphase=3; % number of phases

% phase 1
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=10e3;

% phase 2
GEOM(2).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(2).top=10e3;
GEOM(2).bot=20e3;

% phase 3 (asthenosphere)
GEOM(3).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(3).top=20e3;
GEOM(3).bot=100e3;


% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)

% phase 1: constant viscosity phase
MAT(1).phase=1;
% density parameters
MAT(1).rho0=1000;
MAT(1).alpha=0;
% thermal parameters
MAT(1).k=2;
MAT(1).cp=1000;
% elasticity 
MAT(1).G=1e18;
% diffusion creep parameters
MAT(1).pre_diff=.5/1e18;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
% dislocation creep parameters
MAT(1).pre_disc=.5/1e18;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
% plasticity
MAT(1).mu=0.6;
MAT(1).mumin=0.3;
MAT(1).Cmax=30e6;
MAT(1).Cmin=10e6;
MAT(1).ecrit=0.1;

% phase 2: constant viscosity lithosphere
MAT(2).phase=2;
% density parameters
MAT(2).rho0=3000;
MAT(2).alpha=0;
% thermal parameters
MAT(2).k=2;
MAT(2).cp=1000;
% elasticity 
MAT(2).G=1e18;
% diffusion creep parameters
MAT(2).pre_diff=.5/1e24;
MAT(2).Ediff=0;
MAT(2).ndiff=1;
% dislocation creep parameters
MAT(2).pre_disc=.5/1e24;
MAT(2).Edisc=0;
MAT(2).ndisc=1;
% plasticity
MAT(2).mu=0.6;
MAT(2).mumin=0.3;
MAT(2).Cmax=30e6;
MAT(2).Cmin=10e6;
MAT(2).ecrit=0.1;

% phase 3: constant viscosity astheonsphere
MAT(3).phase=3;
% density parameters
MAT(3).rho0=3000;
MAT(3).alpha=0;
% thermal parameters
MAT(3).k=2;
MAT(3).cp=1000;
% elasticity 
MAT(3).G=1e18;
% diffusion creep parameters
MAT(3).pre_diff=.5/1e20;
MAT(3).Ediff=0;
MAT(3).ndiff=1;
% dislocation creep parameters
MAT(3).pre_disc=.5/1e20;
MAT(3).Edisc=0;
MAT(3).ndisc=1;
% plasticity
MAT(3).mu=0.6;
MAT(3).mumin=0.3;
MAT(3).Cmax=30e6;
MAT(3).Cmin=10e6;
MAT(3).ecrit=0.1;


% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
PARAMS.p0cell=0; % pressure in the top-left corner of the domain (anchor point)

% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant)
% 3/ value of normal velocity 

BC.top=[1 3 0];
BC.bot=[1 0 -2.2182e-10];
BC.left=[1 0 -3.1688e-10];
BC.right=[1 0 3.1688e-10];
PARAMS.BalanceStickyLayer=1; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface


% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 2=Neumann)
% 2/ value
BCtherm.top=[1 0];
BCtherm.bot=[1 1300];
BCtherm.left=[0 0];
BCtherm.right=[0 0];

% TMorrow 22 Sep 2019 - Dike injection routines/controls
BC.DIKE.on=1; % Diking on (1) or off (0)
BC.DIKE.mval=1.0; % M value
BC.DIKE.xL=99.999e3; % left bound of dike
BC.DIKE.xR=100.001e3; % right bound of dike
BC.DIKE.top=10e3; % dike top
BC.DIKE.bot=20e3; % dike bottom
BC.DIKE.injmat=2; % injected phase

% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.YNElast=1; % elasticity on (1) or off (0)
PARAMS.YNPlas=1; % plasticity on (1) or off (0)
PARAMS.tau_heal=1e20; % healing time for plasticity (s)
PARAMS.gx=0; % gravity along x
PARAMS.gy=9.8; % gravity along y
PARAMS.fracCFL=0.5; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=1e25; % maximum viscosity
PARAMS.etamin=1e18; % minimum viscosity
PARAMS.Tsolve=0; % yes (1) or no (0) solve for temperature
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
%PARAMS.a0=-144;
%PARAMS.a1=0.0144;
%PARAMS.a2=0;
%PARAMS.a3=0;
%PARAMS.amp=0; % amplitude of sinusoidal perturbation
%PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=0;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=MAT(1).rho0; 
PARAMS.kref=2;
PARAMS.cpref=1000;

% other thermal parameters
PARAMS.Hradio=0.75e-6; % W/m^3 % radiogenic heat production in radioactive phase
PARAMS.Pradio=[2 4]; % radioactive phase
PARAMS.Nu=10;
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)

zb = (ysize-GEOM(1).bot);
k1 = ((BCtherm.bot(2)-BCtherm.top(2)) + PARAMS.Hradio*zb^2/(2*MAT(1).k))/zb;

PARAMS.a0=PARAMS.T0-k1*GEOM(1).bot-GEOM(1).bot*GEOM(1).bot*PARAMS.Hradio/(2*MAT(1).k);
PARAMS.a1=(k1+GEOM(1).bot*PARAMS.Hradio/MAT(1).k);
PARAMS.a2=-PARAMS.Hradio/(2*MAT(1).k);
PARAMS.a3=0;
PARAMS.amp=100; % amplitude of sinusoidal perturbation
PARAMS.lam=200e3; % wavelength of sinusoidal perturbation
PARAMS.dt_diff=PARAMS.cpref*PARAMS.rhoref/PARAMS.kref/(2/min(GRID.dx)^2 + 2/min(GRID.dy)^2);

% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1000; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
PARAMS.topo_kappa=1e-8; % diffusivity of topography (m^2/s)

% CUSTOM VISCOSITY FUNCTION
%PARAMS.customviscofunction = @(phase,temperature,strain_rate,stress)
%PARAMS.etaMN=18;
%PARAMS.Tc=600;
%PARAMS.Wtan=10;
%PARAMS.etaMX=25;

% Solver iterations
PARAMS.Npicard_min=10; % minimum number of Picard iterations per time step
PARAMS.Npicard_min_interDike=5; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=100; % maximum number of Picard iterations per time step
PARAMS.Npicard_max_interDike=20; % maximum number of inter-Dike Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-3; % convergence criterion
PARAMS.conv_crit_ResL2_interDike=1e-1; % convergence criterion for inter-dike iterations
PARAMS.pitswitch=0; % number of Picard iterations at which the solver switches to quasi-Newton

% SURFACE PROCESSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LEMPARAMS.whichmodel = 0;
% 0 for nothing
% 1 for diffusion (with diffusivity LEMPARAMS.topodiff)
LEMPARAMS.topodiff = 1e-3; % diffusivity of topography (m^2/s)
% 2 for erosion-deposition law of Olive et al. 2014b 
%   (with erosion rate prefactor LEMPARAMS.ED_prefactor 
%    and slope exponent LEMPARAMS.ED_exponent)
LEMPARAMS.ED_exponent = 1.3;
LEMPARAMS.ED_prefactor = 1e-3;
% 3 for stream power incision 
%   (using Topomod routines from J. Braun, B. Kaus, & S. Castelltort)
LEMPARAMS.hsdiff = 1e-9; % diffusivity of hillslopes (m^2/s)
LEMPARAMS.precip = 1; % precipitation rate KEEP SET TO 1 and vary only erodibility
LEMPARAMS.erod = 1e-7; % erodibility coefficient (/yr)
LEMPARAMS.m_param = 0.5; % power law exponent for drainage area
LEMPARAMS.n_param = 1; % power law exponent for slope
LEMPARAMS.distance_in_out_of_page = xsize/3;
LEMPARAMS.npoints_in_out_of_page = 320/2;
LEMPARAMS.npoints_along_x = 960/2;
LEMPARAMS.YNPlot = 0; % plot evolving topography for model #3
LEMPARAMS.YNSed = 1; % yes or no fill sediments up to SedLevel
LEMPARAMS.SedLevel = -5000;
