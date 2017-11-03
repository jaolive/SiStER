function [eta] = SiStER_flow_law_function(type,pre,Ea,n,R,T,epsII,sII,PARAMS)
% [eta] = SiStER_flow_law_function(type,pre,Ea,n,R,T,epsII,sII)
%  gives generic form of ductile creep law 
%  to be used for both diffusion and dislocation creep in
%  get_ductile_rheology functions
%  inputs:
%  type = 'from_stress' or 'from_strain_rate' to set flow law expression
%  with usual parameters: prefactor (pre), activation energy (Ea, J/mol),
%  exponent (n), gas constant (R), temperature (T, in deg C), 
%  second invariant of strain rate (epsII, s^-1) or second invariant
%  of deviatoric stress (sII, Pa)

if strcmp(type,'from_stress')==1
    
    eta = pre.^(-1).*sII.^(1-n).*exp(Ea./(R.*(T+273.15)));
    
    
elseif strcmp(type,'from_strain_rate')==1
    
    eta = pre.^(-1./n).*epsII.^((1-n)./n).*exp((Ea)./(n.*R.*(T+273.15)));
    
    
elseif strcmp(type,'custom')==1
    
    disp('ERROR ? CUSTOM VISCOSITY FUNCTION NOT CURRENTLY DEFINED')
    disp('this feature will be available in a future update.')
    eta = PARAMS.customviscofunction(phase,temperature,strain_rate,stress);
    
else

    disp('ERROR ? flow law is undefined.')

end

