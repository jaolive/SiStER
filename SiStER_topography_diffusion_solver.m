function [topo_new]=SiStER_topography_diffusion_solver(xc,topo_old,dt,K)
% [topo_new]=SiStER_topography_diffusion_solver(xc,topo_old,dt,K)
% very basic explicit diffusion solver to evolve topography

topo_new=topo_old;

dt_surf=0.5*min(diff(xc)).^2/K;
% can't take diff of xc close to edge- if reseeding takes place might be
% too small

if dt_surf<dt
            % DO SEVERAL TEMPERATURE SOLVES UNTIL dt is reached
            nsolve=ceil(dt/dt_surf);
            dt_solve=dt/nsolve;
            topo_solve=topo_old;

            for ksolve=1:nsolve-1
                
               topo_solve(2:end-1)=topo_solve(2:end-1)+dt_solve*K*(1./(0.5*(xc(3:end)+xc(2:end-1)) - 0.5*(xc(2:end-1)+xc(1:end-2))  )).*( (topo_solve(3:end)-topo_solve(2:end-1))./(xc(3:end)-xc(2:end-1))  -  (topo_solve(2:end-1)-topo_solve(1:end-2))./(xc(2:end-1)-xc(1:end-2))  );
                           
            end
           
            disp(['completed ' num2str(nsolve) ' iterations to diffuse topography'])

   
            topo_new=topo_solve;
else
       
    topo_new(2:end-1)=topo_old(2:end-1)+dt*K*(1./(0.5*(xc(3:end)+xc(2:end-1)) - 0.5*(xc(2:end-1)+xc(1:end-2))  )).*( (topo_old(3:end)-topo_old(2:end-1))./(xc(3:end)-xc(2:end-1))  -  (topo_old(2:end-1)-topo_old(1:end-2))./(xc(2:end-1)-xc(1:end-2))  );
                          
end


topo_new(1)=topo_new(2);
topo_new(end)=topo_new(end-1);