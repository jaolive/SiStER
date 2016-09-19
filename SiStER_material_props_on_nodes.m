%=========================================================================
% Get material properties on nodes from advected properties of markers
% G.Ito 8/16 - replaces previous version, which relied on marker viscosity
% for Picard iterations (J.-A.O.)
%=========================================================================

% PHASE PROPORTIONS AT NORMAL AND SHEAR NODES. G.Ito 8/16
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,im);
phase_n=n2interp(1).data;
phase_n=round(phase_n*1e10)/1e10;  %prevents a case in which phase_n>NPhase

[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,im);
phase_s=n2interp(1).data;
phase_s=round(phase_s*1e10)/1e10; %prevents a case in which phase_n>NPhase

% GET MARKER DENSITIES 
[rhom]=SiStER_get_density(im,Tm,MAT);
% pass density to nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom);
rho  = n2interp(1).data;

% GET MARKER ELASTIC PROPERTIES  G.Ito 8/16
[Gm]=SiStER_get_elastic_moduli(im,MAT);
% pass shear modulus to nodes
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,1./Gm);
Gn=1./(n2interp(1).data);
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,1./Gm);
Gs = 1./(n2interp(1).data);   

% PROPERTIES FOR PLASTICITY  G.Ito 8/16
[cohes]=SiStER_get_cohesion(im,ep,MAT); % cohesion depends on plastic strain
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,cohes);
Cohes_n=n2interp(1).data;
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,cohes);
Cohes_s = n2interp(1).data;  

[Mu_s]=SiStER_get_mean_material_prop(cell2mat({MAT([1:Nphase]).mu}),phase_s);
Mu_n=zeros(Ny,Nx);
[Mu_n(2:end,2:end)]=SiStER_get_mean_material_prop(cell2mat({MAT([1:Nphase]).mu}),phase_n(2:end,2:end));

% ADVECTED strain rate invariant G.Ito 8/16
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,epsIIm);
epsII_n=n2interp(1).data;
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,epsIIm);
epsII_s = n2interp(1).data;  


% OLD STRESSES AND PRESSURES G.Ito 8/16
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
sxxOLD(:,:) = n2interp(1).data;
[sxxOLD_s]=SiStER_interp_normal_to_shear_nodes(sxxOLD,dx,dy);

[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,sxym);
sxyOLD(:,:) = n2interp(1).data;  
[sxyOLD_n]=SiStER_interp_shear_to_normal_nodes(sxyOLD);

%MIGHT WANT TO ADVECT PRESSURES (FOR SPEED?) G.Ito 8/16
pold=p; 
ps_old=SiStER_interp_normal_to_shear_nodes(p,dx,dy);

EXYOLD=EXY;
EXXOLD=EXX;
EXX_sOLD=SiStER_interp_normal_to_shear_nodes(EXX,dx,dy);
EXY_nOLD=SiStER_interp_shear_to_normal_nodes(EXY);

%TEMPERATURE ARRAYS NEEDED FOR DUCTILE RHEOLOGY  G.Ito 8/16
Ts=SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
Ts=Ts.data;
Tn=SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,Tm);
Tn=Tn.data;
