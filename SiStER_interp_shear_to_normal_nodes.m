%==============================================================
function [varN]=SiStER_interp_shear_to_normal_nodes(varS)
% 
% G.Ito 8/2016
%===============================================================

varN=zeros(size(varS));
varN(2:end,2:end)=(varS(1:end-1,1:end-1) + varS(2:end,1:end-1)+....
                   varS(1:end-1,2:end)   + varS(2:end,2:end))./4;
return
