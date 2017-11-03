function [eta]=SiStER_get_ductile_rheology_on_nodes_from_strain_rate(MAT,PARAMS,T,epsII,phase_node)
% [eta]=SiStER_get_ductile_rheology_on_nodes_from_strain_rate(MAT,PARAMS,T,epsII,phase_node)
% computes ductile rheology given phase numbers around a node or cell 
% (phase_node is either phase_n or phase_s)
%
% G.Ito 8/2016 % JAO 9/2016 fixed parentheses in exp(Ea/(nRT))
% B.Z. Klein 07/17: can now handle nodes surrounded by any number of phases
% that are NOT consecutively numbered (now looping through all phases)
% thanks to new phase interp functions.

pre_diff = permute(repmat([MAT.pre_diff]', [1,size(T)]), [2,3,1]);
ndiff    = permute(repmat([MAT.ndiff]', [1,size(T)]), [2,3,1]);
Ediff    = permute(repmat([MAT.Ediff]', [1,size(T)]), [2,3,1]);
pre_disc = permute(repmat([MAT.pre_disc]', [1,size(T)]), [2,3,1]);
ndisc    = permute(repmat([MAT.ndisc]', [1,size(T)]), [2,3,1]);
Edisc    = permute(repmat([MAT.Edisc]', [1,size(T)]), [2,3,1]);

epsII = repmat(epsII, [1, 1, PARAMS.Nphase]);
T = repmat(T, [1, 1, PARAMS.Nphase]);


% TO BE REPLACED BY SiStER_flow_law_function SOON
%eta_diff = pre_diff.^(-1./ndiff).*epsII.^((1-ndiff)./ndiff).* ...
%    exp((Ediff)./(ndiff.*PARAMS.R.*(T+273.15)));
%eta_disc = pre_disc.^(-1./ndisc).*epsII.^((1-ndisc)./ndisc).* ...
%    exp((Edisc)./(ndisc.*PARAMS.R.*(T+273.15)));

[eta_diff] = SiStER_flow_law_function('from_strain_rate',pre_diff,Ediff,ndiff,PARAMS.R,T,epsII,zeros(size(epsII)),PARAMS);
[eta_disc] = SiStER_flow_law_function('from_strain_rate',pre_disc,Edisc,ndisc,PARAMS.R,T,epsII,zeros(size(epsII)),PARAMS);


% linearly average between viscosity of each phase type
eta_diffAVG = sum(eta_diff.*phase_node, 3);
eta_discAVG = sum(eta_disc.*phase_node, 3);

eta=(1./eta_diffAVG+1./eta_discAVG).^-1;
eta(eta<PARAMS.etamin)=PARAMS.etamin;

