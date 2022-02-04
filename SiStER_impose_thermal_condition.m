% =========================================================================
% SiStER_impose_thermal_condition
% Imposes one of several thermal conditions on markers, interpolates to
% nodes TMorrow 23 Nov 2021
% =========================================================================

% ===== Thermal Model 1 ===================================================

% Halfspace cooling model set from mean topography, using model spreading
% rate and adding 1 million years of cooling "overlap" in center (ensures
% mantle temperatures do not exist at surface in center of domain)

if BC.DIKE.imposethermal==1    

	% draped on mean topo
	%Tm=BCtherm.bot(2).*erf(0.5*(abs(ym)-mean(topo_y))./sqrt(MAT(2).k/MAT(2).rho0/MAT(2).cp.*(1.0e6*24*3600*365.25+abs(xm-xsize/2)./BC.right(3))));

    % draped on topo
    Tm=BCtherm.bot(2).*erf(0.5*(abs(ym)-interp1(topo_x,topo_y,xm))./sqrt(MAT(2).k/MAT(2).rho0/MAT(2).cp.*(1.0e6*24*3600*365.25+abs(xm-xsize/2)./BC.right(3))));
    
end

% =========================================================================

% ===== Thermal Model 2 ===================================================

% "Notch" structure after Behn & Ito, 2008

if BC.DIKE.imposethermal==2
    
    BC.DIKE.Nwide=40e3;
    BC.DIKE.Lmin=10e3;
    BC.DIKE.Lmax=20e3;
    
    Tm(abs(xm-xsize/2) < BC.DIKE.Nwide)=(ym(abs(xm-xsize/2) < BC.DIKE.Nwide)-mean(topo_y)).*700.*((1/BC.DIKE.Lmax-1/BC.DIKE.Lmin).*abs(xm(abs(xm-xsize/2) < BC.DIKE.Nwide)-xsize/2)./BC.DIKE.Nwide+1/BC.DIKE.Lmin);
	Tm(abs(xm-xsize/2) >= BC.DIKE.Nwide)=(ym(abs(xm-xsize/2) >= BC.DIKE.Nwide)-mean(topo_y)).*700.*((1/BC.DIKE.Lmax-1/BC.DIKE.Lmin)+1/BC.DIKE.Lmin);
	
end

% =========================================================================

% bound temperature
Tm(Tm>BCtherm.bot(2))=BCtherm.bot(2);
Tm(Tm<0)=0;

% air is T=0
Tm(im<=1.5)=0;

% get temperature on nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
T(:,:)  = n2interp(1).data;
