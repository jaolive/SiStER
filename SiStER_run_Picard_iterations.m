% SiStER_run_Picard_iterations

for pit=1:PARAMS.Npicard_max

    disp(['t = ' num2str(t) ': Picard iteration ' num2str(pit) ' out of ' num2str(PARAMS.Npicard_max)]); %G.Ito
    
    SiStER_flow_solve;  
    % ASSESS CONVERGENCE:
    % convergence is considered achieved when 
        
    % get strain rate from current solution on normal nodes
    EXYc=(EXY(1:Ny-1,1:Nx-1)+EXY(2:Ny,1:Nx-1)+EXY(1:Ny-1,2:Nx)+EXY(2:Ny,2:Nx))/4;
    EPSII=sqrt(EXX(2:Ny,2:Nx).^2+EXYc.^2); 
    depsII=abs(EPSII(:)-EPSIIold(:))./(max(EPSIIold(:))); % relative change in strain rate
    % fraction of domain size where change in strain rate is greater than
    % PARAMS.conv_crit1:
    frac_depsII=sum(depsII>PARAMS.conv_crit1)/(length(depsII)); 
   
    % get velocity field from current solution on shear nodes
    SiStER_interp_velocities_to_shear_nodes;  
    v=sqrt(vxc.^2+vyc.^2); % total velocity at shear nodes
    dv=abs(v(:)-vold(:))./median(vold(:)); % relative change in velocity
    % fraction of domain size where change in velocity field is greater than
    % PARAMS.conv_crit1:
  	frac_dv=sum(dv>PARAMS.conv_crit1)/(Nx*Ny);   
    
    % average fraction of large changes
    frac_change=.5*(frac_dv+frac_depsII);
    
    % convergence assessment
    if (frac_change<= PARAMS.conv_crit2 && pit >= PARAMS.Npicard_min)  %both strain rate and volume below crit1 at local minima
        disp([num2str(pit) ' Picard iterations converged: relative solution change is greater than ' num2str(100*PARAMS.conv_crit1) '% in only ' num2str(100*frac_change) ' % of the domain.']); % JAO
        break;
	elseif (pit==PARAMS.Npicard_max);
       	disp([num2str(pit) ' Picard iterations failed to converge.']);
    end

    % Save solutions from this Picard iteration to compare to the next.
	vold=v;   
	EPSIIold=EPSII;

end

      
