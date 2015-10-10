%function [indexes] = nn(scanA, scanB, thres)
%Compute the nearst neighbour of each point in scan A from scan B
%
%INPUTS:
% scanA a N x 2 matrix with scan values for A 
% scanB a N x 2 matrix with scan values for B
% thres a threshold on maximum distance
%
%OUTPUTS:
% a matrix Nmatches x 2 with the corresponce indexes of the two scans
%
% Implemented by Fabio Tozeto Ramos 18/08/05

function [indtemp] = nn(scanA, scanB, thres)
warning off
if nargin<3, thres = 0.5; end
[NA] = size(scanA);
[NB] = size(scanB);
%indexes=zeros(2,361);
XA2 = sum(scanA.^2,2);
XB2 = sum(scanB.^2,2);
XB2T = XB2';
distance = XA2(:,ones(1,NB))+XB2T(ones(1,NA),:)-2*scanA(:,1)*scanB(:,1)'-2*scanA(:,2)*scanB(:,2)';
[mindist,indtemp] = min(distance,[],2);
indtemp = [(1:NA)' indtemp];
ind = (mindist<thres^2);
indtemp = indtemp(ind,:);
