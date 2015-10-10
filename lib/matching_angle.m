%function [indexes] = matching_angle(scanA,scanB,thres)
%Compute the nearst neighbour of each point in scan A from scan B
%
%INPUTS:
% scanA a N x 2 matrix with scan values for A 
% scanB a N x 2 matrix with scan values for B
% thres a threshold on maximum angle difference
%
%OUTPUTS:
% a matrix Nmatches x 2 with the corresponce indexes of the two scans
%
% Implemented by Fabio Tozeto Ramos 18/08/05
function [indexes] = matching_angle(scanA,scanB,thres)

if nargin<3, thres = deg2rad(15); end % 15 degrees

[NA,D] = size(scanA);
[NB,D] = size(scanB);
%Computes the polar coordenates

% Angle
angleA = atan(scanA(:,2)./scanA(:,1));
angleB = atan(scanB(:,2)./scanB(:,1));
%Fix -2*pi, +2*pi problems
%angleA = pi_to_pi(angleA);
%angleB = pi_to_pi(angleB);
angleBT = angleB';

% Radius
radiusA = sqrt(sum(scanA.^2,2));
radiusB = sqrt(sum(scanB.^2,2));
radiusBT = radiusB';

% Computes the module of the difference between the radius of a point
% in radiusA and in radiusB.
diff_radius = sqrt((radiusA(:,ones(1,NB))-radiusBT(ones(1,NA),:)).^2);
[mindist_r,ind_r] = min(diff_radius');
diff_angle = sqrt((angleA-angleB(ind_r)).^2);
indexes=find(diff_angle<abs(thres));
indexes=[indexes ind_r(indexes)'];
%[mindist,indexes] = min(diff_angle');
%indexes = [[1:NA]' indexes'];
%ind = find(mindist<thres);
%indexes = indexes(ind,:);
