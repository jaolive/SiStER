function [X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid(xsize,ysize,GRID)
% [X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid(xsize,ysize,GRID)
% sets up the staggered grid coordinates

%% Coded for both regular and variable grid spacing G.Ito 7/15
%
%GRID.x defines the points in x within which grid spacing can differ
%for x<=GRID.x(1), the APPROXIMATE grid spacing is GRID.dx(1); 
%for GRID.x(1) < x <= GRID.x(2), the approx. spacing is GRID.dx(2);
%for GRID.x(end-1) < x <= xsize, the approx. spacing is GRID.dx(end);
%
%Thus GRID.dx MUST be longer than GRID.x by 1 entry.
%GRID.x=[0] and GRID.dx=[0.5 1]*1e3 will sent a UNIFORM spacing of 1000m
%The grid is defined in SiStER_initialize_grid_GI.  This is also 
%where Nx will be defined.
%
%To ensure an EXACT grid spacing, be sure GRID.dx fits as an integer number
%of times in each grid interval GRID.x

%% The same logic is used for y.


x=0;
y=0;
for i=1:length(GRID.x);
  n=round((GRID.x(i)-x(end))/GRID.dx(i));  %set number of element in zone i
  dd=(GRID.x(i)-x(end))/n;                 %set real grid size
  x=[x, x(end)+(1:n).*dd];
end
n=round((xsize-x(end))./GRID.dx(end));
dd=(xsize-x(end))./n;
x=[x, x(end)+(1:n).*dd];

for i=1:length(GRID.y);
  n=round((GRID.y(i)-y(end))./GRID.dy(i));
  dd=(GRID.y(i)-y(end))./n;
  y=[y, y(end)+(1:n).*dd];
end
n=round((ysize-y(end))./GRID.dy(end));
dd=(ysize-y(end))./n;
y=[y, y(end)+(1:n).*dd];

Nx=length(x);
Ny=length(y);  % they may have changed !

% 2D Grid
[X,Y] = meshgrid(x,y);

% create vectors of grid spacing
dx=diff(x);
dy=diff(y);

% Create vectors for cell centers positions (staggered nodes)
xc=0.5*(x(2:Nx)+x(1:Nx-1));
yc=0.5*(y(2:Ny)+y(1:Ny-1));

disp(['** Nx = ' num2str(Nx) ' Ny = ' num2str(Ny) ' ***']);


