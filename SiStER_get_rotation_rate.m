function [om]=SiStER_get_rotation_rate(vx,vy,dx,dy,BC)
% [om]=SiStER_get_rotation_rate(vx,vy,dx,dy,BC)
% computes the rotation rate (vorticity) at shear nodes


[Ny, Nx]=size(vx);
om=zeros(size(vy));


for i=1:Ny
    for j=1:Nx   
        
        if i>=2 && i<=Ny-1 && j>=2 && j<=Nx-1
            om(i,j)=0.5*(-2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
        end
        
        if i==1 && j>=2 && j<=Nx-1 % top
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end         
            om(i,j)= 0.5*(-dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif i==Ny && j>=2 && j<=Nx-1 % bottom
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end         
            om(i,j)= 0.5*(-dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif j==1 && i>=2 && i<=Ny-1
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end     
            om(i,j)=0.5*(-2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        elseif j==Nx && i>=2 && i<=Ny-1
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end     
            om(i,j)=0.5*(-2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        end
        
        % corners
        if i==1 && j==1
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end   
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end   
            om(i,j)=(-dvxdy + dvydx)/2;
            
        elseif i==1 && j==Ny
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end  
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end  
            om(i,j)=(-dvxdy + dvydx)/2;
            
        elseif i==Nx && j==Ny
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end  
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end  
            om(i,j)=(-dvxdy + dvydx)/2;
            
        elseif i==Nx && j==1
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end  
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end  
            om(i,j)=(-dvxdy + dvydx)/2;
            
        end
            
            
        
    end
end        