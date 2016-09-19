function [eta]=SiStER_get_ductile_rheology(MAT,PARAMS,T,epsII,phase_node)
% 
% computes ductile rheology given phase numbers around a node or cell 
% (phase_node is either phase_n or phase_s)
% Potential shortcoming:  it assumes only <+ TWO phases that are consecutively
% numbered surround each node
% G.Ito 8/2016 % JAO 9/2016 fixed parentheses in exp(Ea/(nRT))
%-------------------------------------------------------------------------
         
[Ny,Nx]=size(epsII);

ph1=floor(phase_node);
ph2=ph1;

if (sum(ceil(phase_node)>ph1+1));
    disp('problems in SiStER_get_ductile_rheology');
    halt
end
pre_diff=cell2mat({MAT(1:end).pre_diff})'; 
ndiff=cell2mat({MAT(1:end).ndiff})';
Ediff=cell2mat({MAT(1:end).Ediff})';
pre_disc=cell2mat({MAT(1:end).pre_disc})'; 
ndisc=cell2mat({MAT(1:end).ndisc})';
Edisc=cell2mat({MAT(1:end).Edisc})';

eta_diff1=pre_diff(ph1).^(-1./ndiff(ph1)).*epsII.^((1-ndiff(ph1))./ndiff(ph1)).*...
          exp(Ediff(ph1)./(ndiff(ph1).*PARAMS.R.*(T+273.15)));
eta_disc1=pre_disc(ph1).^(-1./ndisc(ph1)).*epsII.^((1-ndisc(ph1))./ndisc(ph1)).*...
          exp(Edisc(ph1)./(ndisc(ph1).*PARAMS.R.*(T+273.15)));

eta_diff2=zeros(Ny,Nx);
eta_disc2=zeros(Ny,Nx);
f1=ones(Ny,Nx);

ii=find(phase_node > ph1);
if (length(ii)>1);
ph2(ii)=ph1(ii)+1;
eta_diff2(ii)=pre_diff(ph2(ii)).^(-1./ndiff(ph2(ii))).*epsII(ii).^((1-ndiff(ph2(ii)))./ndiff(ph2(ii))).*...
    exp(Ediff(ph2(ii))./(ndiff(ph2(ii)).*PARAMS.R.*(T(ii)+273.15)));
eta_disc2(ii)=pre_disc(ph2(ii)).^(-1./ndisc(ph2(ii))).*epsII(ii).^((1-ndisc(ph2(ii)))./ndisc(ph2(ii))).*...
    exp(Edisc(ph2(ii))./(ndisc(ph2(ii)).*PARAMS.R.*(T(ii)+273.15)));
f1(ii)=(ph2(ii)-phase_node(ii));  %fraction of phase 1
end;


eta=f1.*(1./eta_diff1+1./eta_disc1).^-1+(1-f1).*(1./eta_diff2+1./eta_disc2).^-1;
    
% if (0);  %testing the shape of cell2mat
%     YY=[[1:10]' [11:20]'];
%     stuff(1).g=1;
%     stuff(2).g=2;
%     ph=YY;
%     ph(:,1)=1; ph(:,2)=2;
%     check1=cell2mat({stuff(ph).g});
%     ZZ=YY(:)'.*check1;
%     ZZ=reshape(ZZ,10,2);
% end
return
