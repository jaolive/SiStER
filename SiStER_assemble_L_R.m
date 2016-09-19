function [L, R, Kc, Kb]=SiStER_assemble_L_R(dx,dy,etas,etan,rho,BC,PARAMS,srhs_xx,srhs_xy)


%% Fill LHS and RHS Solution Matrices

p0cell=PARAMS.p0cell;

Nx=length(dx)+1;
Ny=length(dy)+1;

gx=PARAMS.gx;
gy=PARAMS.gy;

% Vector of right part initialization
R=zeros(Nx*Ny*3,1);
Lii = zeros(10*Nx*Ny*3,1);
Ljj = zeros(10*Nx*Ny*3,1);
Lvv = zeros(10*Nx*Ny*3,1);
nn = 1;

%%%%%%%%%% SCALING THE FD MATRIX
% Computing Kc and Kb coefficients
meta=min(min(etas));
mdx=max(dx);
mdy=max(dy);
Kc=2*meta/(mdx+mdy);
Kb=4*meta/(mdx+mdy)^2;


% pressure anchor
%IP=2;
%JP=3;
JP=ceil((Nx-1)/2); %G.Ito
IP=2;


% fill in FD matrix and right-hand side

for j=1:Nx
    for i=1:Ny
    %for j=1:Nx
        
        in=(j-1)*Ny+i;
        inp=3*in-2;
        invx=3*in-1;
        invy=3*in;
        
        
        %%%%%%% CONTINUITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (i==1) || (j==1) || (i==2 && j==2) || (i==2 && j==Nx) || ....
           (i==Ny && j==2) || (i==Ny && j==Nx) || (i==IP && j==JP && BC.top(2)~=3) || ...
           (BC.top(2)==3 && i==2 && j > 2 && j < Ny) %<--G.Ito
            
            
            % boundary conditions
            
            if (i==1) || (j==1) % P(i,j)=0
                %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            if (i==2 && j==2) || (i==Ny && j==2)  % left corners P(i,j)-P(i,j+1)=0
                %                L(inp,inp)=Kb;
                %                L(inp,inp+3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp+3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            if (i==2 && j==Nx) || (i==Ny && j==Nx)  % right corners P(i,j)-P(i,j-1)=0
                %                L(inp,inp)=Kb;
                %                L(inp,inp-3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp-3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            
            if (BC.top(2)==3 && i==2 && j > 2 && j < Ny)  %open top effect all nodes but corners, G.Ito
                % Pressure gradient between top two rows extrapolates to 0 pressure at very top
% 
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb.*(1+dy(i-1)/(dy(i)+dy(i-1)));
                nn = nn+1;

                Lii(nn) = inp;
                Ljj(nn) = inp+3;
                Lvv(nn)  = -Kb.*dy(i-1)/(dy(i)+dy(i-1));
                nn = nn+1;

                R(inp,1)=0;
               
            end
            
            if (i==IP && j==JP && BC.top(2)~=3)   % additional pressure cell i=2, j=3 P(i,j)=p0cell
                 %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(inp,1)=Kb*p0cell/Kc;
                           
            end
           
            
        else
            
            
            % internal nodes
            
            % coeffs for vx
            %             L(inp,invx-3)=Kc/dx(j-1);
            %             L(inp,invx-3*Ny-3)=-Kc/dx(j-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3;
            Lvv(nn)  = Kc/dx(j-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3*Ny-3;
            Lvv(nn)  = -Kc/dx(j-1);
            nn = nn+1;
            
            % coeffs for vy
            %             L(inp,invy-3*Ny)=Kc/dy(i-1);
            %             L(inp,invy-3*Ny-3)=-Kc/dy(i-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = Kc/dy(i-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny-3;
            Lvv(nn)  = -Kc/dy(i-1);
            nn = nn+1;
            
            
            % RHS
            R(inp,1)=0;%-drhodt(i,j);
            
        end
        
        
        %%%%%%% X-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if ((j==1) && (i<=Ny-1)) || ((j==Nx) && (i<=Ny-1)) || (i==1 && j<=Nx-1 && j>=2) || (i==Ny-1 && j<=Nx-1 && j>=2) || (i==Ny)
            
            
            
            
            if ((j==1) && (i<=Ny-1)) %|| ((j==Nx) && (i<=Ny-1))  % X-STOKES:  left boundary 
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.left(3);
            end
            
            if ((j==Nx) && (i<=Ny-1))  %  X-STOKES:  right  boundary
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.right(3);
            end
            
            
            if i==1 && j<=Nx-1 && j>=2 % X-STOKES:  upper boundary not including sides
                if BC.top(1)==1 % free slip
                    %                  L(invx,invx+3)=Kc/(0.5*(dy(i)+dy(i+1)));
                    %                  L(invx,invx)=-Kc/(0.5*(dy(i)+dy(i+1)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = -Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                elseif BC.top(1)==0  % no slip
                    %                  L(invx,invx+3)=Kc/(dy(i)+dy(i+1));
                    %                  L(invx,invx)=-Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    % Multiply whole eq. by 2 so coefficients are closer to
                    % the others
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = -2*Kc/(dy(i)+dy(i+1));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.top_profile(j);  %G.Ito x-velocity at top boundary;
                end
                
            end
            
            if i==Ny-1 && j<=Nx-1 && j>=2 % X-STOKES: lower boundary not including sides
                if BC.bot(1)==0 % no slip
                    %                  L(invx,invx)=Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(dy(i-1)+dy(i));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.bot_profile(j);
                    
                elseif BC.bot(1)==1 % free slip
                    %                  L(invx,invx)=Kc/(0.5*(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(0.5*(dy(i-1)+dy(i)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                end
                
                
                
            end
            
            if i==Ny  % lower boundary, ghost vx=0
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=0;
            end
            
            
            
        else

            % X-STOKES: internal nodes
            
            % X-STOKES: coeffs for vx
            
            %         L(invx,invx-3*Ny)=4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3*Ny;
            Lvv(nn)  = 4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx)=(-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx;
            Lvv(nn)  = (-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3*Ny)=4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx-3)=2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3)=2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1)));
            nn = nn+1;
            
            % X-STOKES: coeffs for vy
            
            %         L(invx,invy+3)=2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy+3-3*Ny)=-2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3-3*Ny;
            Lvv(nn)  = -2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy)=-2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy;
            Lvv(nn)  = -2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy-3*Ny)=2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            % X-STOKES: coeffs for Pp
            
            %         L(invx,inp+3*Ny+3)=-2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            %         L(invx,inp+3)=2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3;
            Lvv(nn)  = 2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            % RHS
            R(invx,1)=-gx*(rho(i,j)+rho(i+1,j))/2 - ...
                2.*(srhs_xx(i+1,j+1)-srhs_xx(i+1,j))/(dx(j)+dx(j-1)) - ...
                (srhs_xy(i+1,j)-srhs_xy(i,j))/dy(i);
            
        end
        
        
        
        
        %%%%%%% Y-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (j==1 && i<=Ny-1 && i>=2) || (j==Nx-1 && i<=Ny-1 && i>=2) || (j==Nx) || (i==1 && j<=Nx-1) || (i==Ny && j<=Nx-1)
                   
            if j==1 && i>=2 && i<=Ny-1  % left boundary without edges
                if BC.left(1)==1 % free slip
                    %              L(invy,invy)=-2*Kc/(dx(j)+dx(j+1));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j)+dx(j+1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.left(1)==0 % no slip
                    %              L(invy,invy)=-2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j+1)+dx(j));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j+1)+dx(j));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            if j==Nx-1 && i<=Ny-1 && i>=2  % right boundary without edges 
                if BC.right(1)==1 % free slip
                    %              L(invy,invy)=Kc/(0.5*(dx(j)+dx(j-1)));
                    %              L(invy,invy-3*Ny)=-Kc/(0.5*(dx(j)+dx(j-1)));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn)  = -Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.right(1)==0 % no slip
                    %              L(invy,invy)=Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    %              L(invy,invy-3*Ny)=-Kc/(dx(j)+dx(j-1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn) = -Kc/(dx(j)+dx(j-1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
                
            end
            if j==Nx  % right boundary, ghost vy=0
                %              L(invy,invy)=Kb;
                
                Lii(nn) = invy;
                Ljj(nn) = invy;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invy,1)=0;
                
            end
            
            if i==1 && j<=Nx-1  % upper boundary
                if BC.top(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.top(3);
                elseif BC.top(2)==3
                    %               open top dvy/dy=0, G.Ito
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kb;
                    nn      = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3;
                    Lvv(nn) = -Kb;
                    nn      = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            
            if i==Ny && j<=Nx-1  % lower boundary without edges
                
                if BC.bot(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.bot(3);
                    
                elseif BC.bot(2) == 2 %%  external boundary
                    
                    dL=300e3;
                                        
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kc*(1/dL+1/dy(i-1));
                    nn = nn+1;
                    
                    Lii(nn) = invy; 
                    Ljj(nn) = invy-3;
                    Lvv(nn) = Kc*(-1/dy(i-1));
                    nn = nn+1;

                    R(invy, 1) = 0;
                end
                    
            end
            
            
        else
            
            
            
            % internal nodes
            
            % coeffs for vy
            %         L(invy,invy-3*Ny)=2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy)=(-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy;
            Lvv(nn)  = (-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy+3*Ny)=2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1)));
            nn = nn+1;
            
            %         L(invy,invy-3)=4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3;
            Lvv(nn)  = 4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i)));
            nn  = nn+1;
            
            %         L(invy,invy+3)=4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            % coeffs for vx
            
            %         L(invy,invx+3*Ny)=2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));  %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx+3*Ny-3)=-2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny-3;
            Lvv(nn)  = -2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            
            %         L(invy,invx)=-2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx;
            Lvv(nn)  = -2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx-3)=2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn=nn+1;
            
            % coeffs for Pp
            %         L(invy,inp+3*Ny+3)=-2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            %         L(invy,inp+3*Ny)=2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny;
            Lvv(nn)  = 2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            % RHS (note: sxx = -syy because deviatoric)
            R(invy,1)=-gy*(rho(i,j)+rho(i,j+1))/2 + ...
                2.*(srhs_xx(i+1,j+1)-srhs_xx(i,j+1))/(dy(i)+dy(i-1)) - ...
                (srhs_xy(i,j+1)-srhs_xy(i,j))/dx(j); %
            
        end
        

        
    end % for j
end % for i

nn = nn-1;
%%% Build Sparse Matrix

Li = Lii(1:nn);
Lj = Ljj(1:nn);
Lv = Lvv(1:nn);


L = sparse(Li,Lj,Lv);
% SOLVE LINEAR SYSTEM now down in SiStER_flow_solve G.Ito
% Matlab direct solver
% tic
% S=L\R;
% toc
% 





end

