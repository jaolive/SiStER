%==========================================================================
% SiStER_flow_solve
% Performs inner solve of linear LS=R system as well as outer, iterative
% solution for non-linear dependence of L (viscosity) on S (vx,vy,P)
% Used to be named "run_Picard_iterations" but name changed by G.Ito 6/21/16
%==========================================================================


if PARAMS.BalanceStickyLayer==1
% BALANCE FLUXES %%% JAO July 16, 2015
% RE-ADJUST BCs SO FLUX OF ROCK AND STICKY AIR MATERIAL ARE CONSERVED
% locate height of sticky layer - rock contact on the sides
    bL=interp1(topo_x,topo_y,0);
    bR=interp1(topo_x,topo_y,xsize);
    utop=BC.right(3)*(bL+bR)/xsize;
    ubot=BC.right(3)*(2*ysize-bL-bR)/xsize;
    BC.top(3)=utop;
    BC.bot(3)=-ubot;
end


ResL2=1; 

for pit=1:PARAMS.Npicard_max
    
    
    if pit==1
        ResL2init=ResL2;
    end

    %% ---------------------------------------------------------------------------------
    % Compute visco-elasto-plastic viscosities
    %---------------------------------------------------------------------------------
    SiStER_VEP_rheology   

    %---------------------------------------------------------------------------------
    % Assemble L and R matrices
    %---------------------------------------------------------------------------------
    [L, R, Kc, Kb]=SiStER_assemble_L_R(dx,dy,Zs.*etas,Zn.*etan,rho,BC,PARAMS,srhs_xx,srhs_xy); %G.Ito
    
    %---------------------------------------------------------------------------------
    % Residual:  L and R are from current solution S
    %---------------------------------------------------------------------------------
    if (exist('S','var'));
        Res=L*S-R;
        ResL2=norm(Res,2)/norm(R,2);
        
    else
        
        S=L\R; disp('First solve is linearized');
        Res=L*S-R;
        ResL2=norm(Res,2)/norm(R,2);
        
        
    end
    %---------------------------------------------------------------------------------
    % Solve for new solution S using Picard or approximate Newton or a
    % combination of the two 
    %---------------------------------------------------------------------------------  
    
    if (pit >= PARAMS.pitswitch);
       if pit==PARAMS.pitswitch; disp('switching from Picard to approx. Newton'); end;
       beta=1;
       S=S-beta.*(L\Res);  % approximate Newton update, with L as approximation to Jacobian
       it_type='Newton: ';
    else
       S=L\R; % Picard update
       it_type='Picard: ';
    end
   
    [p, vx, vy]=SiStER_reshape_solver_output(S,Kc,Nx,Ny);
        
    %% ASSESS CONVERGENCE    
    if(ResL2<PARAMS.conv_crit_ResL2 && pit >= PARAMS.Npicard_min)
        disp(['Final residual = ' num2str(ResL2)])
        disp([num2str(pit) ' iterations converged: L2 norm of residual dropped below ' num2str( PARAMS.conv_crit_ResL2)]);
        break;
	elseif (pit==PARAMS.Npicard_max);
        disp(['Final residual = ' num2str(ResL2)])
        disp(['WARNING! ' num2str(pit) ' Picard / approx. Newton iterations failed to converge within tolerance of ' num2str( PARAMS.conv_crit_ResL2)]);
    end

%end
   
%% get strain rate on nodes current solution
[EXX,EXY]=SiStER_get_strain_rate(vx,vy,dx,dy,BC);

EXY_n=SiStER_interp_shear_to_normal_nodes(EXY);       
EXX_s=SiStER_interp_normal_to_shear_nodes(EXX,dx,dy); 
epsII_n=sqrt(EXX.^2+EXY_n.^2);
epsII_s=sqrt(EXX_s.^2+EXY.^2);



% helpful to visualize convergence
% figure(1)
% pcolor(X,Y,log10(etas))
% set(gca,'ydir','reverse')
% axis equal
% caxis([18 25])
% colorbar
% title(num2str(pit))
% pause(.001)

% RESIDUAL FOR INDIVIDUAL VARIABLES
% [pres, vxres, vyres]=SiStER_reshape_solver_output(Res,Kc,Nx,Ny);
% figure(1)
% pcolor(X,Y,vxres)
% set(gca,'ydir','reverse')
% axis equal
% colorbar
% title(num2str(pit))
% pause(.001)


end

      
