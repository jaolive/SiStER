% BISECTION ALGORITHM FOR LOCAL ITERATIONS
% calculates a deviatoric stress sII that can be fed to the creep laws
% and is consistent with the current viscosity
% this appears to be a necessary step to use a stress-based creep law
% when elasticity is on (i.e., the strain rate is not entirely viscous)

% STARTING VISCOSITY (TO BE ADJUSTED)
visco=etas_new;

% INITIALIZE
sIILOW=sIIOLD_s;
sIIHIGH=sIIOLD_s;
[visco0]=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Ts,sIIOLD_s,phase_s);
if PARAMS.YNElast==1
    Zs=Gs.*dt_m./(Gs.*dt_m+visco0);
else
    Zs=ones(Ny,Nx);
end
sIINEW=2*visco0.*epsII_s.*Zs+(1-Zs).*sIIOLD_s;
sIIHIGH(sIINEW>sIIHIGH)=sIINEW(sIINEW>sIIHIGH);
sIILOW(sIINEW<sIILOW)=sIINEW(sIINEW<sIILOW);

[viscoNEW]=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Ts,sIINEW,phase_s);
if PARAMS.YNElast==1
    Zs=Gs.*dt_m./(Gs.*dt_m+viscoNEW);
else
    Zs=ones(Ny,Nx);
end 
sIINEW=2*viscoNEW.*epsII_s.*Zs+(1-Zs).*sIIOLD_s;
sIIHIGH(sIINEW>sIIHIGH)=sIINEW(sIINEW>sIIHIGH);
sIILOW(sIINEW<sIILOW)=sIINEW(sIINEW<sIILOW);

% 10 iterations of a BISECTION ALGORITHM
Nbisec=10;
for i=1:Nbisec
    
    sII=0.5*(sIIHIGH+sIILOW);
    sIIPAST=sII;
    [visco]=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,Ts,sII,phase_s);
    if PARAMS.YNElast==1
        Zs=Gs.*dt_m./(Gs.*dt_m+visco);
    else
        Zs=ones(Ny,Nx);
    end
    sII=2*visco.*epsII_s.*Zs+(1-Zs).*sIIOLD_s;

    sIILOW(sII>sIIPAST)=sIIPAST(sII>sIIPAST);
    sIIHIGH(sII<=sIIPAST)=sIIPAST(sII<=sIIPAST);

end
