function [EXX,EXY]=SiStER_get_strain_rate(vx,vy,dx,dy,BC)
% [EXX EXY]=SiStER_get_strain_rate(vx,vy,dx,dy,BC)
% computes the deviatoric strain rate: EXX, at normal nodes
% and EXY, at shear nodes, from vx and vy on the eulerian grid
%
% G.Ito 11/12/15 updated calculations at corner nodes to be
%    consistent with the top and bottom boundary conditions, as those
%    (not the side) conditions control the solutions for vx,vy at the
%    points closest to the corners in SiStER_Stokes_solver.m

[Ny, Nx]=size(vx);

EXX=zeros(size(vx));
EXY=zeros(size(vy));

% Getting EXX (on normal nodes) is straightforward
for i=2:Ny
    for j=2:Nx   
        EXX(i,j)=(vx(i-1,j)-vx(i-1,j-1))/dx(j-1);        
    end
end
        
% Getting EXY (on shear nodes) requires some thinking when dealing with
% boundaries
for i=1:Ny
    for j=1:Nx   
        
        if i>=2 && i<=Ny-1 && j>=2 && j<=Nx-1
            EXY(i,j)=0.5*(2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
        end
        
        if i==1 && j>=2 && j<=Nx-1 % top
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end         
            EXY(i,j)= 0.5*(dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif i==Ny && j>=2 && j<=Nx-1 % bottom
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                %dvxdy=-2*vx(i-1,j)/dy(i-1);
                dvxdy=2.*(BC.bot(4)-vx(i-1,j))/dy(i-1);  %G.Ito needed for tangential velocity on base
            end         
            EXY(i,j)= 0.5*(dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif j==1 && i>=2 && i<=Ny-1 
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end     
            EXY(i,j)=0.5*(2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        elseif j==Nx && i>=2 && i<=Ny-1
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end     
            EXY(i,j)=0.5*(2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        end
        
        % corners
        if i==1 && j==1
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end   
            
            if BC.top(2)==0 % vy imposed on top (Dirichlet)
                dvydx=0;
            elseif BC.top(2)==3 % open top so set to be the same as adjacent node
                dvydx=2*(vy(i,j+1)-vy(i,j))/(dx(j)+dx(j+1));
            else
                disp(['Stopped in get_strain_rate:  EXY for corner elements not coded for BC.top(2)='....
                num2str(BC.top(2))]);
                halt
            end   
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        elseif i==1 && j==Nx 
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no slip
                dvxdy=2*vx(i,j)/dy(i);
            end  
            if BC.top(2)==0 % vy imposed on top (Dirichlet)
                dvydx=0;
            elseif BC.top(2)==3 % open top so set to be the same as adjacent node
                 dvydx=2*(vy(i,j-1)-vy(i,j-2))./(dx(j-1)+dx(j-2));
            end  
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        elseif i==Ny && j==Nx
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end
            if BC.bot(2)==0 % vy imposed on bottom (Dirichlet)
                dvydx=0;
            elseif BC.bot(2)==3 % open bottom so set to be the same as adjacent node
                dvydx=2*(vy(i,j-1)-vy(i,j-2))./(dx(j-1)+dx(j-2));
            else
                disp(['Stopped in get_strain_rate:  EXY for corner elements not coded for BC.bot(2)='....
                num2str(BC.bot(2))]);
                halt
            end
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        elseif i==Ny && j==1
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end  
            if BC.bot(2)==0 % vy imposed on bottom (Dirichlet)% rollers
                dvydx=0;
            elseif BC.bot(2)==3 % open bottom so set to be the same as adjacent node
                dvydx=2*(vy(i,j+1)-vy(i,j))./(dx(j+1)+dx(j));
            end  
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        end
            
            
        
    end
end        