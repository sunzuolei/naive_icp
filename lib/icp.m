%function [pos,indexes,err]=icp(rawscanA, rawscanB,  relpos, thres_trans, nit, interp, plotflag)
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

function [pos,indexes,err]=icp(rawscanA, rawscanB,  relpos, thres_trans, nit, interp,plotflag)

if nargin<5; thres_trans = 0.5; interp=0; nit=50; end
if nargin<6; interp=0; end

if size(rawscanA,1)~=2
    scanA = get_laser_points(rawscanA, 70); % [x; y]; (2*361)
    scanB = get_laser_points(rawscanB, 70);
else
    scanA = rawscanA;
    scanB = rawscanB;
end
pos = relpos;
prevpos=relpos;
for i=1:nit
    if ~interp % if don't interpolation
        scanAglobal = transform_to_global(scanA, pos); % scanAglobal (2*361)
        indexes=nn_mutual(scanAglobal',scanB',thres_trans);
        % a matrix Nmatches x 2 with the corresponce indexes of the two scans
        %     indexes(:,1) - the index of point in ScanA which can find its correspondence in ScanB
        %     indexes(:,2) - the index of point in ScanB which correspond the
        %     points indicated by indexes(:,1)
        [post(1),post(2),post(3)]=rel_pose(scanAglobal(:,indexes(:,1))',scanB(:,indexes(:,2))');
        pos=pos+post';
    else % if interpolation
        scanAglobal = transform_to_global(scanA, pos);
        [indexes, intscan]=nn_interp(scanAglobal',scanB',thres_trans);
        [post(1),post(2),post(3)]=rel_pose(scanAglobal(:,indexes(:,1))',intscan);
        pos=pos+post';
    end
    
    err = sum(sum((scanAglobal(:,indexes(:,1))-scanB(:,indexes(:,2))).^2));
    % err = sum[(XA1-XB1)^2+(YA1-YB1)^2,...,(XAn-XBn)^2+(YAn-YBn)^2];
    %Stop if converged
    posdist=sqrt((pos(1)-prevpos(1))^2+(pos(2)-prevpos(2))^2);
    if posdist<1e-8
        break;
    else
        prevpos=pos;
    end
end
%scA=scanA(:,indexes(:,1));
%scB=scanB(:,indexes(:,2));

% figure;hold on;
% plot(scanA(1,:),scanA(2,:),'b+');
% plot(scanB(1,:), scanB(2,:),'r+');
% 
if plotflag
    figure;hold on;
    scanAglobal =  transform_to_global(scanA, pos);
    plot(scanAglobal(1,:), scanAglobal(2,:),'b+');
    plot(scanB(1,:),scanB(2,:),'r+');
end
% Transforms laser ranges to cartesians coordenates
% by TIm Bailey
function x = get_laser_points(r, MAXR)
phi = -pi/2:(pi/360):pi/2;
ii = find(r < MAXR);
r = r(ii); phi = phi(ii);
x = [r.*cos(phi); r.*sin(phi)];
