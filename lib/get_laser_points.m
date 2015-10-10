%function x = get_laser_points(r, MAXR)
%Transform the laser's range into x,y coordenates
function x = get_laser_points(r, MAXR)
phi = -pi/2:(pi/360):pi/2;
ii = find(r < MAXR);
r = r(ii); phi = phi(ii);
x = [r.*cos(phi); r.*sin(phi)];
