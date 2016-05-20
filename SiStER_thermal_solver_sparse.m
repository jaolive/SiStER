function [T, rhs, Lii, Ljj, Lvv]=SiStER_thermal_solver_sparse(x,y,Told,rho,cp,kx,ky,dt,BCtherm)
% [T, rhs, Lii, Ljj, Lvv]=SiStER_thermal_solver_sparse(x,y,Told,rho,cp,kx,ky,dt,BCtherm)
% implicit Solver for thermal diffusion
% rho cp dT/dt = div (k grad T)
% J.-A. Olive, November 2014
% B.Z. Klein, added sparse matrix filling, End of 2014


% need to use sparse-filling for L

Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

Li=zeros(Nx*Ny, 1);
Lj=zeros(Nx*Ny, 1);
Lv=zeros(Nx*Ny, 1);
rhs=zeros(Nx*Ny,1);

n = 1;

for i=1:Ny
    for j=1:Nx
        
        in=i+Ny*(j-1); % global index
        
        if i==1 % top boundary
            
            if BCtherm.top(1)==1 % Dirichlet
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n)  = 1;
                
                n = n+1;
               
                rhs(in) = BCtherm.top(2);
                
            %L(in,in)=1;
            %rhs(in)=BCtherm.top(2);
            
            elseif BCtherm.top(1)==2 % Neumann
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = -1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in+1;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in)=BCtherm.top(2)*dy(i);
                
%             L(in,in)=-1;
%             L(in,in+1)=1;
%             rhs(in)=BCtherm.top(2)*dy(i);    
            
            end
            
        elseif i==Ny % bottom boundary
            
            
            if BCtherm.bot(1)==1 % Dirichlet
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in) = BCtherm.bot(2);
                
%             L(in,in)=1;
%             rhs(in)=BCtherm.bot(2);
            
            elseif BCtherm.bot(1)==2 % Neumann
                
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in-1;
                Lv(n) = -1;
                
                n = n+1;
                
                rhs(in) = BCtherm.bot(2)*dy(i-1);
                
%             L(in,in)=1;
%             L(in,in-1)=-1;
%             rhs(in)=BCtherm.bot(2)*dy(i-1);    
            
            end
 
        elseif j==1 % left boundary
            
            
            if BCtherm.left(1)==2 % Neumann
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = -1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in+Ny;
                Lv(n) = 1;
                
                n = n+1;
            
                rhs(in) = BCtherm.left(2)*dx(j);
                
%             L(in,in)=-1;
%             L(in,in+Ny)=1;
%             rhs(in)=BCtherm.left(2)*dx(j);  
            
            elseif BCtherm.left(1)==1 % Dirichlet
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in) = BCtherm.left(2);
                
%             L(in,in)=1;
%             rhs(in)=BCtherm.left(2); 
            
            end
            
            
        elseif j==Nx % right boundary
            
            
            if BCtherm.right(1)==2 % Neumann
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in-Ny;
                Lv(n) = -1;
                
                n = n+1;
                
                rhs(in) = BCtherm.right(2)*dx(j-1);
                
%             L(in,in)=1;
%             L(in,in-Ny)=-1;
%             rhs(in)=BCtherm.right(2)*dx(j-1);  
            
            elseif BCtherm.right(1)==1 % Dirichlet
                
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in) = BCtherm.right(2);
                
%             L(in,in)=1;
%             rhs(in)=BCtherm.right(2); 
            
            end
            
            
            
        else
            
            
%         %   internal nodes
            
            ddx=dx(j-1)+dx(j);
            ddy=dy(i-1)+dy(i);
            
            Li(n) = in;
            Lj(n) = in;
            Lv(n) = rho(i,j)*cp(i,j)+2*dt*kx(i,j)/(dx(j)*ddx) + 2*dt*kx(i,j-1)/(dx(j-1)*ddx) + ...
                2*dt*ky(i,j)/(ddy*dy(i)) + 2*dt*ky(i-1,j)/(ddy*dy(i-1));
            
            n = n+1;
            
%             L(in,in)=rho(i,j)*cp(i,j)+2*dt*kx(i,j)/(dx(j)*ddx) + 2*dt*kx(i,j-1)/(dx(j-1)*ddx) + ...
%                 2*dt*ky(i,j)/(ddy*dy(i)) + 2*dt*ky(i-1,j)/(ddy*dy(i-1));
            
            Li(n) = in;
            Lj(n) = in+Ny;
            Lv(n) = -2*dt*kx(i,j)/(dx(j)*ddx);
            
            n = n+1;
            
%             L(in,in+Ny)=-2*dt*kx(i,j)/(dx(j)*ddx);
% 
            Li(n) = in;
            Lj(n) = in-Ny;
            Lv(n) = -2*dt*kx(i,j-1)/(dx(j-1)*ddx);
            
            n = n+1;
% 
% %             L(in,in-Ny)=-2*dt*kx(i,j-1)/(dx(j-1)*ddx);
%             
            Li(n) = in;
            Lj(n) = in+1;
            Lv(n) = -2*dt*ky(i,j)/(dy(i)*ddy);
            
            n = n+1;
% 
% %             L(in,in+1)=-2*dt*ky(i,j)/(dy(i)*ddy);
% 
            Li(n) = in;
            Lj(n) = in-1;
            Lv(n) = -2*dt*ky(i-1,j)/(dy(i-1)*ddy);
            
            n = n+1;

%             L(in,in-1)=-2*dt*ky(i-1,j)/(dy(i-1)*ddy);
            
            rhs(in)=rho(i,j)*cp(i,j)*Told(i,j);
        
        end
        
    end
end

nn = n-1;

Lii = Li(1:nn);
Ljj = Lj(1:nn);
Lvv = Lv(1:nn);

L = sparse(Lii, Ljj, Lvv);

tvec=L\rhs;

T=reshape(tvec,Ny,Nx);

    

    
        
        
        
        
