function [etam,syield]=SiStER_get_viscosity(im,Tm,epsII,MAT,PARAMS,ep,press)
% [etam,syield]=SiStER_get_viscosity(im,Tm,epsII,MAT,PARAMS,ep,press)
% Computes the viscosity and yield stress as a function of material type
% (im), temperature (Tm, in degrees C), second invariant of strain rate 
% (epsII), accumulated plastic strain (ep), pressure (press), 
% and various parameters stored in MAT & PARAMS.
%
% DEFAULT CREEP LAW
% eta for diffusion and dislocation creep have the form:
% eta = pre^(-1/n) * epsII^((1-n)/n) * exp(E/nRT);
% then both get harmonically averaged, and harmonically averaged 
% with a plastic viscosity, if needed
%
% J.-A. Olive, summer 2015 (first cut of current version)
% Accelerated by B. Klein & G.Ito 11/15 with a loop on marker types 
% (using "unique" function)


type=unique(im);

eta_diff=zeros(size(im));
eta_disc=zeros(size(im));

for i=1:length(type);
    MATprops = MAT(type(i));
    logic = im == type(i);
    eta_diff(logic)=[MATprops.pre_diff].^(-1./[MATprops.ndiff]).*epsII(logic).^((1-[MATprops.ndiff])./[MATprops.ndiff]).*exp([MATprops.Ediff]./([MATprops.ndiff].*PARAMS.R.*(Tm(logic)+273.15)));
    eta_disc(logic)=[MATprops.pre_disc].^(-1./[MATprops.ndisc]).*epsII(logic).^((1-[MATprops.ndisc])./[MATprops.ndisc]).*exp([MATprops.Edisc]./([MATprops.ndisc].*PARAMS.R.*(Tm(logic)+273.15)));
end;


[syield] = SiStER_get_yield_stress(im,ep,press,MAT);
eta_plas = syield./(2*epsII);


if PARAMS.YNPlas==0
    etam=(1./eta_diff+1./eta_disc + 1/PARAMS.etamax).^-1;
elseif PARAMS.YNPlas==1
    etam=(1./eta_diff+1./eta_disc + 1/PARAMS.etamax + 1./eta_plas).^-1;
end

% enforce lower-limit on viscosity
etam(etam<PARAMS.etamin)=PARAMS.etamin; 

