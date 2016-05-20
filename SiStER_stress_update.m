% update elastic stresses on markers following a solve (but before advection)

% get marker rotation rates
[ROT]=SiStER_get_rotation_rate(vx,vy,dx,dy,BC);
[om]=SiStER_interp_shear_nodes_to_markers(ROT,x,y,xm,ym,icn,jcn);

% compute stress changes on nodes
dsxx=(2*etan.*EXX-sxxOLD).*Zn;
dsxy=(2*etas.*EXY-sxyOLD).*Zs;

% pass stress changes to markers
[dsxxm]=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,xm,ym,icn,jcn);
[dsxym]=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,xm,ym,icn,jcn);

% UPDATE STRESSED ON MARKERS
sxxm=sxxm+dsxxm;
sxym=sxym+dsxym;

% ROTATE MARKERS
alpha=om*dt_m;
sxymtemp = sxxm.*(sin(alpha).^2) + sxym.*cos(2.*alpha);
sxxm = sxxm.*(cos(alpha).^2 - sin(alpha).^2)- sxym.*sin(2.*alpha);
sxym=sxymtemp;

% get the current second invariant of strain rate to update ep below
% (not needed if PARAMS.epsII_from_stress==0, because we already have it 
% from the last Picard iteration)
if (PARAMS.YNPlas==1 && PARAMS.epsII_from_stress==1)
    [exxm]=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn);
    [exym]=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn);
    epsIIm=sqrt(exxm.^2+exym.^2);
end

% PLASTICITY UPDATE
if PARAMS.YNPlas==1
   % STRESS REDUCTION IN YIELDED REGIONS
   sIItemp=sqrt(sxxm.^2+sxym.^2);
   sxxm(sIItemp>=.99*syield)=sxxm(sIItemp>=.99*syield).*syield(sIItemp>=.99*syield)./sIItemp(sIItemp>=.99*syield); 
   sxym(sIItemp>=.99*syield)=sxym(sIItemp>=.99*syield).*syield(sIItemp>=.99*syield)./sIItemp(sIItemp>=.99*syield); 
   % STRAIN ACCUMULATION IN YIELDED REGIONS
   ep(sIItemp>=.99*syield)=ep(sIItemp>=.99*syield)+dt_m*epsIIm(sIItemp>=.99*syield);
   epNH(sIItemp>=.99*syield)=ep(sIItemp>=.99*syield)+dt_m*epsIIm(sIItemp>=.99*syield); % non-healed plastic strain
   % STRAIN HEALING
   ep=ep./(dt_m/PARAMS.tau_heal+1);
end



