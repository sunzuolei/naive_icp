%function pos=icp_car(scanA, scanB,  relpos, nit, thres_trans)
%Implements the ICP algorithm with no interpolation
%
%INPUTS:
% rawscanA, rawscanB vectors with laser scans
% relpos intial relative position
% thres_trans threshold on maximum translation
% nit number of iterations
% interp flag to perform interpolation
%
%OUTPUTS:
% pos relative position
%
% Implemented by Fabio Tozeto Ramos

function [pos,scA ,scB,err]=icp_car(scanA, scanB,  relpos, thres_trans, nit, interp)

if nargin<5, thres_trans = 0.5, interp=0, nit=50; end

%scanA = get_laser_points(rawscanA, 40);
%scanB = get_laser_points(rawscanB, 40);
pos = relpos;
prevpos=relpos;
for i=1:nit
    if ~interp
        scanAglobal = transform_to_global(scanA, pos);
        indexes=nn(scanAglobal',scanB',thres_trans);
        [post(1),post(2),post(3)]=rel_pose(scanAglobal(:,indexes(:,1))',scanB(:,indexes(:,2))');
        pos=pos+post';
    else
        scanAglobal = transform_to_global(scanA, pos);
        [indexes, intscan]=nn_interp(scanAglobal',scanB',thres_trans);
        [post(1),post(2),post(3)]=rel_pose(scanAglobal(:,indexes(:,1))',intscan);
        pos=pos+post';
    end
    
    err = sum(sum((scanAglobal(:,indexes(:,1))-scanB(:,indexes(:,2))).^2));
    %Stop if converged
    posdist=sqrt((pos(1)-prevpos(1))^2+(pos(2)-prevpos(2))^2);
    if posdist<1e-8
        break;
    else
        prevpos=pos;
    end
end
scA=scanA(:,indexes(:,1));
scB=scanB(:,indexes(:,2));

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
