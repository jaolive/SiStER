%=========================================================================
% Computes velocities at cell corners (shear nodes) from vx and vy 
% on staggered grid 
% (it can be useful to have vx and vy at the same location)
% G.Ito 7/15
%=========================================================================


%% ------------------------------------------------------------------------
% Vx is defined on the x values of cell edges so only need to interpolate
% to y values of cell tops and bottoms
% -------------------------------------------------------------------------
vxc=zeros(Ny,Nx);
% internal nodes
vxc(2:Ny-1,:)=vx(1:Ny-2,:).*(((1-dy(1:Ny-2)./(dy(1:Ny-2)+dy(2:Ny-1))))'*ones(1,Nx))+....
              vx(2:Ny-1,:).*(((1-dy(2:Ny-1)./(dy(1:Ny-2)+dy(2:Ny-1))))'*ones(1,Nx));

% Top 
if (BC.top(1)==1);  %free slip
    vxc(1,:)=vx(2,:);
elseif (BC.top(1)==0); %No slip
    vxc(1,:)=0;
else
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.top(1) ~= 0 or 1');
    halt
end;
% Bottom
if (BC.bot(1)==1);  %free slip
    vxc(Ny,:)=vx(Ny-1,:);
    
elseif (BC.bot(1)==0); %No slip
    if (length(BC.bot)==4);
      vxc(Ny,:)=BC.bot(4);
    else
      vxc(Ny,:)=0;
    end;
else
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.bot(1) ~= 0 or 1');
    halt
end;



%% ------------------------------------------------------------------------
% Vy is defined on y values of cell tops and bottoms so only need to 
% interpolate to x-values of cell edges
% -------------------------------------------------------------------------
vyc=zeros(Ny,Nx);
% internal values
vyc(:,2:Nx-1)=vy(:,1:Nx-2).*(ones(Ny,1)*(1-dx(1:Nx-2)./(dx(1:Nx-2)+dx(2:Nx-1))))+....
              vy(:,2:Nx-1).*(ones(Ny,1)*(1-dx(2:Nx-1)./(dx(1:Nx-2)+dx(2:Nx-1))));

%Left
if (BC.left(1)==1); %free slip
    vyc(:,1)=vy(:,2);
elseif (BC.left(1)==0); % no slip
    vyc(:,1)=0;
else
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.left(1) ~= 0 or 1');
    halt
end;

%Right
if (BC.right(1)==1); %free slip
    vyc(:,Nx)=vy(:,Nx-1);
    
elseif (BC.right(1)==0); % no slip
    vyc(:,Nx)=0.0;
    
else 
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.right(1) ~= 0 or 1');
    halt
end; 
    


%

