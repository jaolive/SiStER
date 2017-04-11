function [quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy)
% [quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy)
% Tells a marker which cell (and which quadrant of that cell) it belongs to.
    % icn,jcn are the indexes of the upper-left shear node of the cell
    % that a given marker is currently in
    % quad is the quadrant of the cell that contains the marker
    % (quad = 1 means bottom-right, then numbered clockwise)
    % sped up by B. Klein in Fall 2016 by using interp1 function


%% Determine Location of Markers and Quadrant of Element 
M=length(xm);
icn = zeros(1,M);
jcn = icn;
quad = icn; % quadrant 1 = bottom-right, numbered clockwise

indX = 1:length(x);
indY = 1:length(y);

Ix = interp1(x, indX, xm, 'nearest', 'extrap');
Iy = interp1(y, indY, ym, 'nearest', 'extrap');

% [~, Ix] = min(abs(bsxfun(@minus, xm, x')));
% [~, Iy] = min(abs(bsxfun(@minus, ym, y')));

jcn(xm>x(Ix))  = Ix(xm>x(Ix));
jcn(xm<=x(Ix)) = Ix(xm<=x(Ix))-1;

icn(ym>y(Iy))  = Iy(ym>y(Iy));
icn(ym<=y(Iy)) = Iy(ym<=y(Iy))-1;

jcn(jcn==0) = 1;
jcn(jcn>length(dx)) = length(dx);

icn(icn==0) = 1;
icn(icn>length(dy)) = length(dy);

disx = abs((xm-x(jcn))./dx(jcn));
disy = abs((ym-y(icn))./dy(icn));


xRIGHT = disx > 0.5;
yUP = disy > 0.5;

quad(xRIGHT & yUP) = 3;
quad(xRIGHT & ~yUP)  = 2;
quad(~xRIGHT & yUP)  = 4;
quad(~xRIGHT & ~yUP)   = 1;
