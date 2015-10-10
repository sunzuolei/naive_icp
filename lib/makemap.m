%Produce a map from laser and odometry data
%makemap(LSR1,DR)
function makemap(LSR1,DR)

figure, axis equal, axis([-350 500 -450 450]); grid, hold on
laser=LSR1.Range;
% Configurable values
INTERP = 0; % 0 - simple ICP, 1 - ICP with interpolation
NIT = 50;   % number of ICP iterations for each scan-match
GATE1 = 3; % 1st data association gate for each point in scan
GATE2 = 0.5;% 2nd data association gate for each point in scan

% Constants
MAXR = 70; % this laser data times-out at ~7.9 m

% Initialise
x = [0;0;0];
delta = [0;0;0];
%p1 = get_laser_points(laser(1,2:end), MAXR);
p1=laser(1,:);
pgp=get_laser_points(p1, MAXR);
% Loop
for i=2:size(laser,1)
    %p2 = get_laser_points(laser(i,2:end), MAXR);
    p2 = laser(i,:); disp(i);
    %delta = scan_matching_LM(p1, p2, delta, GATE1, 15)
    delta = odometrymov(DR,[0;0;0],LSR1.Time(i-1),LSR1.Time(i));
    delta= -delta;
    delta = icp(p1, p2, delta, GATE1, 15,0,0);
    [delta, indexes] = icp(p1, p2, delta, GATE2, NIT, INTERP,0);
    %delta = idc(p1, p2, delta, GATE2, NIT);
    
    x = transform_to_global(-delta, x);
    
    % plots
    pt = get_laser_points(p2, MAXR);
    pg = transform_to_global(pt, x);
    plot(x(1), x(2), '+', pg(1,:), pg(2,:), '.','erasemode','none');
    drawnow
    pt = get_laser_points(p1, MAXR);
    plot_assoc_gen(pgp,pg,indexes(:,1),indexes(:,2));
    p1 = p2;
    pgp = pg;
end

%
%

function x = get_laser_points(r, MAXR)
phi = -pi/2:(pi/360):pi/2;
ii = find(r < MAXR);
r = r(ii); phi = phi(ii);
x = [r.*cos(phi); r.*sin(phi)];

function plot_assoc_gen(scA, scB, Apts, Bpts)
pa = make_association_lines(scA, scB, Apts, Bpts);
pt = make_association_lines(scA, scB, Apts, Bpts);
%plot(scA(1,:),scA(2,:),'r+', scB(1,:),scB(2,:),'b+');
if ~isempty(pa)
    plot(pa(1,:),pa(2,:),'g');
end
if ~isempty(pt)
    plot(pt(1,:),pt(2,:),'m');
end
function p = make_association_lines(a, b, Apts,Bpts)
p = [];
for i=1:length(Apts)
    if Bpts(i)~=0
        p = [p a(:, Apts(i)) b(:,Bpts(i)) [NaN;NaN]];
    end
end