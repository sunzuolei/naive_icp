%function [indexes] = nnT(scanA, scanB, thres)
%Compute the nearst neighbour of each point in scan A from scan B
%
%INPUTS:
% scanA a 2 x N matrix with scan values for A 
% scanB a 2 x N matrix with scan values for B
% thres a threshold on maximum distance
%
%OUTPUTS:
% a matrix Nmatches x 2 with the corresponce indexes of the two scans
%
% Implemented by Fabio Tozeto Ramos 18/08/05

function [indtemp] = nnT(scanA, scanB, thres)

if nargin<3, thres = 0.5; end
[D,NA] = size(scanA);
[D,NB] = size(scanB);
indexes=zeros(2,351);
XA2 = sum(scanA.^2);
XB2 = sum(scanB.^2);
XB2T = XB2';
distance = XA2(ones(1,NB),:)+XB2T(:,ones(1,NA))-2*scanB(1,:)'*scanA(1,:)-2*scanB(2,:)'*scanA(2,:);
[mindist,indtemp] = min(distance');
indtemp = [[1:NA]' indtemp'];
ind = find(mindist<thres^2);
indtemp = indtemp(ind,:);
