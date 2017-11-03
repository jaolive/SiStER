function [eta]=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,T,sII,phase_node)
% [eta]=SiStER_get_ductile_rheology_on_nodes_from_stresses(MAT,PARAMS,T,sII,phase_node)
% computes ductile rheology given phase numbers around a node or cell 
% (phase_node is either phase_n or phase_s)
%
% G.Ito 8/2016 % JAO 9/2016 fixed parentheses in exp(Ea/(nRT))
% B.Z. Klein 07/17: can now handle nodes surrounded by any number of phases
% that are NOT consecutively numbered (now looping through all phases)
% thanks to new phase interp functions.
% written by J.A. Olive based on BZ Klein's
% get_ductile_rheology_on_nodes_from_strain_rate.m function (Nov. 2017)
%-------------------------------------------------------------------------

pre_diff = permute(repmat([MAT.pre_diff]', [1,size(T)]), [2,3,1]);
ndiff    = permute(repmat([MAT.ndiff]', [1,size(T)]), [2,3,1]);
Ediff    = permute(repmat([MAT.Ediff]', [1,size(T)]), [2,3,1]);
pre_disc = permute(repmat([MAT.pre_disc]', [1,size(T)]), [2,3,1]);
ndisc    = permute(repmat([MAT.ndisc]', [1,size(T)]), [2,3,1]);
Edisc    = permute(repmat([MAT.Edisc]', [1,size(T)]), [2,3,1]);

sII(sII==0)=1e-3; % to avoid zero stress at first time step

sII = repmat(sII, [1, 1, PARAMS.Nphase]);
T = repmat(T, [1, 1, PARAMS.Nphase]);

% TO BE REPLACED BY SiStER_flow_law_function SOON
%eta_diff=pre_diff.^(-1).*sII.^(1-ndiff).*...
%         exp(Ediff./(PARAMS.R.*(T+273.15)));
%eta_disc=pre_disc.^(-1).*sII.^(1-ndisc).*...
%          exp(Edisc./(PARAMS.R.*(T+273.15)));
           
[eta_diff] = SiStER_flow_law_function('from_stress',pre_diff,Ediff,ndiff,PARAMS.R,T,zeros(size(sII)),sII,PARAMS);
[eta_disc] = SiStER_flow_law_function('from_stress',pre_disc,Edisc,ndisc,PARAMS.R,T,zeros(size(sII)),sII,PARAMS);

            
% linearly average between viscosity of each phase type
eta_diffAVG = sum(eta_diff.*phase_node, 3);
eta_discAVG = sum(eta_disc.*phase_node, 3);

eta=(1./eta_diffAVG+1./eta_discAVG).^-1;
eta(eta<PARAMS.etamin)=PARAMS.etamin;

return
