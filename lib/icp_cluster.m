%function pos=icp_cluster(rawscanA, rawscanB,  relpos, nit, thres_trans)
%Implements the ICP algorithm over clusters
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
% Implemented by Fabio Tozeto Ramos 09/10/2006

function pos=icp_cluster(rawscanA, rawscanB,  relpos, thres_trans, nit)

if nargin<5, thres_trans = 0.5, nit=50; end

scanA = get_laser_points(rawscanA, 60);
scanB = get_laser_points(rawscanB, 60);
pos = relpos;
%Cluster points
[clusterA]=cluster_laserV2(scanA,1,0);
[clusterB]=cluster_laserV2(scanB,1,0);
maxA=max(clusterA);
maxB=max(clusterB);
meansA=zeros(2,maxA);
meansB=zeros(2,maxB);
maxclust=max([max(clusterA) max(clusterB)]);
%Find means for each cluster
for j=1:maxclust
    if j<=maxA
        meansA(:,j)=mean(scanA(:,find(clusterA==j))');
    end
    if j<=maxB
        meansB(:,j)=mean(scanB(:,find(clusterB==j))');
    end
end
prevpos=relpos;

%Start icp over clusters' means
for i=1:nit
    meansAglobal = transform_to_global(meansA, pos);
    indexes=nn(meansAglobal',meansB',thres_trans);
    [pos(1),pos(2),pos(3)]=rel_pose(meansAglobal(:,indexes(:,1))',meansB(:,indexes(:,2))');
    %Stop if converged
    posdist=sqrt((pos(1)-prevpos(1))^2+(pos(2)-prevpos(2))^2);
    if posdist<0.01
        break;
    else
        prevpos=pos;
    end
end

%figure;hold on;
%plot(scanA(1,:),scanA(2,:),'b+');
%plot(scanB(1,:), scanB(2,:),'r+');

%figure;hold on;
%scanAglobal =  transform_to_global(scanA, pos);
%plot(scanAglobal(1,:), scanAglobal(2,:),'b+');
%plot(scanB(1,:),scanB(2,:),'r+');


% Transforms laser ranges to cartesians coordenates
% by TIm Bailey
function x = get_laser_points(r, MAXR)
phi = -pi/2:(pi/360):pi/2;
ii = find(r < MAXR);
r = r(ii); phi = phi(ii);
x = [r.*cos(phi); r.*sin(phi)];
