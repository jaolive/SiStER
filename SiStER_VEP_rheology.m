%==========================================================================
% SiStER_UPDATE_RHEOLOGY
% calculates visco-elasto-plastic rheology terms on shear and normal nodes
% G.Ito 8/20
%==========================================================================

%--------------------------------------------------------------------------
% plastic viscosity first (Mohr-Coulomb)
%--------------------------------------------------------------------------
if (PARAMS.YNPlas==1);
    dp=p-pold;
    ps=ps_old+SiStER_interp_normal_to_shear_nodes(dp,dx,dy);  %<<<< Interpolation done here
    ps(ps<0)=0;
    pnn=p; 
    pnn(p<0)=0;
    yield_s=(Cohes_s+Mu_s.*ps).*cos(atan(Mu_s));
    yield_n=(Cohes_n+Mu_n.*pnn).*cos(atan(Mu_n));
    eta_plas_s=0.5.*yield_s./max(epsII_s-(yield_s-sqrt(sxxOLD_s.^2+sxyOLD.^2))./(2.*Gs.*dt_m),min(epsII_s(:))*1e-6);  
	eta_plas_n=0.5.*yield_n./max(epsII_n-(yield_n-sqrt(sxxOLD.^2+sxyOLD_n.^2))./(2.*Gn.*dt_m),min(epsII_n(:))*1e-6); 
end

%-------------------------------------------------------------------------
% Ductile viscosity
%-------------------------------------------------------------------------
[etas_new]=SiStER_get_ductile_rheology(MAT,PARAMS,Ts,epsII_s,phase_s);
[etan_new(2:end,2:end)]=SiStER_get_ductile_rheology(MAT,PARAMS,Tn(2:end,2:end),epsII_n(2:end,2:end),phase_n(2:end,2:end));


if PARAMS.YNPlas==1 
    

    % identify yielding nodes to update ep
    s_nodes_yield=find(etas_new>eta_plas_s);
    n_nodes_yield=find(etan_new>eta_plas_n);
    % incorporate plastic viscosity into effective viscosity
    etan_new(2:end,2:end)=(1./eta_plas_n(2:end,2:end) + 1./etan_new(2:end,2:end) + 1./PARAMS.etamax).^-1;
    etas_new=(1./eta_plas_s + 1./etas_new + 1./PARAMS.etamax).^-1;

    
    % possible alternative = sharp viscosity drop where yielding occurs (may be less stable)
    % etan_new(n_nodes_yield)=eta_plas_n(n_nodes_yield);
    % etas_new(s_nodes_yield)=eta_plas_s(s_nodes_yield);

else
    
    etan_new(2:end,2:end)=(1./etan_new(2:end,2:end) + 1./PARAMS.etamax).^-1;
    etas_new=(1./etas_new + 1./PARAMS.etamax).^-1;

end


etan=etan_new;
etas=etas_new;
etan(etan<PARAMS.etamin)=PARAMS.etamin;
etas(etas<PARAMS.etamin)=PARAMS.etamin;


%-------------------------------------------------------------------------
% ELASTICITY TERMS 
%-------------------------------------------------------------------------
if PARAMS.YNElast==0
    Zs=ones(size(etas));
    Zn=ones(size(etan));
else
    Zs=Gs*dt_m./(Gs*dt_m+etas);
    Zn=Gn*dt_m./(Gn*dt_m+etan);
end
% right-hand size (1-Z)*sigmaOLD_ij
srhs_xx=(1-Zn).*sxxOLD;
srhs_xy=(1-Zs).*sxyOLD;

