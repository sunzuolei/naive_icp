% Implements the Iterative Dual Correspondence Algorithm (IDC)

function pos=idc(rawscanA, rawscanB,  relpos, thres_trans, nit)

if nargin<5, thres_trans = 0.5; end
thres_rot = deg2rad(25);
scanA = get_laser_points(rawscanA, 40);
scanB = get_laser_points(rawscanB, 40);
pos = relpos;
prevpos=relpos;
post = zeros(3,1);
posr = zeros(3,1);
for i=1:nit
    scanAglobal = transform_to_global(scanA, pos);
    indexes=nn(scanAglobal',scanB',thres_trans);
    [post(1),post(2),post(3)]=rel_pose(scanAglobal(:,indexes(:,1))',scanB(:,indexes(:,2))');
    scanAglobal = transform_to_global(scanA, post);
    indexes=matching_angle(scanAglobal',scanB',thres_rot);
    if size(indexes,1)>5
        [posr(1),posr(2),posr(3)]=rel_pose(scanAglobal(:,indexes(:,1))',scanB(:,indexes(:,2))');
    end
    pos(1) =pos(1) + posr(1) + post(1);
    pos(2) = pos(2) + posr(2) + post(2);
    pos(3) = pos(3) + posr(3) + post(3);
    
    err = sum(sum((scanAglobal(:,indexes(:,1))-scanB(:,indexes(:,2))).^2))
    %Stop if converged
    posdist=sqrt((pos(1)-prevpos(1))^2+(pos(2)-prevpos(2))^2);
    if posdist<1e-8
        break;
    else
        prevpos=pos;
    end
    %pos
    %thres_rot = pos(3);
end
figure;hold on;
plot(scanA(1,:),scanA(2,:),'b+');
plot(scanB(1,:), scanB(2,:),'r+');

figure;hold on;
scanAglobal =  transform_to_global(scanA, pos);
plot(scanAglobal(1,:), scanAglobal(2,:),'b+');
plot(scanB(1,:),scanB(2,:),'r+');

% Transforms laser ranges to cartesians coordenates
% by TIm Bailey
function x = get_laser_points(r, MAXR)
phi = -pi/2:(pi/360):pi/2;
ii = find(r < MAXR);
r = r(ii); phi = phi(ii);
x = [r.*cos(phi); r.*sin(phi)];
