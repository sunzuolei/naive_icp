%function pos=scan_matching_LM(rawscanA, rawscanB,  relpos, nit, thres_trans)
%Implements scan-matching with Levenberg-Marquardt optimisation 
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

function [pos, indexes]=scan_matching_LM(rawscanA, rawscanB,  relpos, thres_trans, nit)

if nargin<5, thres_trans = 0.02, interp=0, nit=50; end
global indexes;
scanA = get_laser_points(rawscanA, 60);
scanB = get_laser_points(rawscanB, 60);
pos = relpos;
prevpos=relpos;

%Cluster points
% [clusterA]=cluster_laserV2(scanA,1,0);
% [clusterB]=cluster_laserV2(scanB,1,0);
% maxA=max(clusterA);
% maxB=max(clusterB);
% meansA=zeros(2,maxA);
% meansB=zeros(2,maxB);
% maxclust=max([max(clusterA) max(clusterB)]);
% %Find means for each cluster
% for j=1:maxclust
%     if j<=maxA
%         meansA(:,j)=mean(scanA(:,find(clusterA==j))');
%     end
%     if j<=maxB
%         meansB(:,j)=mean(scanB(:,find(clusterB==j))');
%     end
% end


warning off;
%Optimisation with the Levenberg Marquardt method
options = optimset('LargeScale','on');
%options = optimset(options, 'LevenbergMarquardt','on');
options = optimset(options, 'MaxFunEvals', 10000);
options = optimset(options, 'MaxIter', 500);
options = optimset(options, 'Display', 'iter');
options = optimset(options, 'TolFun', 1e-3);
options = optimset(options, 'TolX', 1e-8);
options = optimset(options, 'Diagnostics', 'on');
[pos,resnorm,resi,exit] = lsqnonlin(@(pos)scan_match_err(pos, scanA, scanB, thres_trans), relpos,[-2;-2;-pi/4],[2;2;pi/4],options);
%[pos,resnorm,resi,exit] = fminsearch(@(pos)scan_match_err(pos, scanA, scanB, thres_trans), relpos,[],[],options)

% figure;hold on;
% plot(scanA(1,:),scanA(2,:),'b+');
% plot(scanB(1,:), scanB(2,:),'r+');
% 
% figure;hold on;
% scanAglobal =  transform_to_global(scanA, pos);
% plot(scanAglobal(1,:), scanAglobal(2,:),'b+');
% plot(scanB(1,:),scanB(2,:),'r+');

% Error function for LM optimisation
function [err,indexes] = scan_match_err(pos, scanA, scanB, thres_trans)
global indexes;
scanAglobal = transform_to_global(scanA, pos);
indexes=nn(scanAglobal',scanB',thres_trans);
err = sum((scanAglobal(:,indexes(:,1))-scanB(:,indexes(:,2))).^2,2);

% Error function for LM optimisation
function err = scan_match_err2(pos, scanA, scanB, thres_trans)
global indexes;
L=2; %free parameter
scanAglobal = transform_to_global(scanA, pos);
indexes=nn(scanAglobal',scanB',thres_trans);
dx=scanB(1,indexes(:,2))-scanAglobal(1,indexes(:,1));
dy=scanB(2,indexes(:,2))-scanAglobal(2,indexes(:,1));
phi=sum(((dx.*scanAglobal(2,indexes(:,1))-dy.*scanAglobal(1,indexes(:,1))).^2)...
    ./(scanAglobal(2,indexes(:,1)).^2+scanAglobal(1,indexes(:,1)).^2+L^2));
err = [sum((scanAglobal(:,indexes(:,1))-scanB(:,indexes(:,2))).^2,2); -phi];

% Error function for LM optimisation
function [err,indexes] = scan_match_err3(pos, scanA, scanB, thres_trans)
global indexes;
scanAglobal = transform_to_global(scanA, pos);
indexes=nn_mutual(scanAglobal',scanB',thres_trans);

err = sum((scanAglobal(:,indexes(:,1))-scanB(:,indexes(:,2))).^2,2);

% Transforms laser ranges to cartesians coordenates
% by TIm Bailey
function x = get_laser_points(r, MAXR)
phi = -pi/2:(pi/360):pi/2;
ii = find(r < MAXR);
r = r(ii); phi = phi(ii);
x = [r.*cos(phi); r.*sin(phi)];
