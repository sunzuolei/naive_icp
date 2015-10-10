%function [indexes] = nn_interp(scanA, scanB, thres)
%Compute the nearst neighbour of each point in scan A from scan B 
%doing interpolation.
%
%INPUTS:
% scanA a N x 2 matrix with scan values for A 
% scanB a N x 2 matrix with scan values for B
% thres a threshold on maximum distance
%
%OUTPUTS:
% a matrix Nmatches x 2 with the corresponce indexes of the two scans
%     indexes(:,1) - index of point in scan A
%     indexes(:,2) - index of point in scan B
%
% Implemented by Fabio Tozeto Ramos 07/10/06

function [indexes, intscaB] = nn_interp(scanA, scanB, thres)

if nargin<3, thres = 0.5; end

[NA,D] = size(scanA); % NA = N
[NB,D] = size(scanB);

%Computes the line-point distance in a vectorised manner
distX=scanB(2:NB,1)'-scanB(1:(NB-1),1)'; % consider using diff. the size is [1*(N-1)]
distY=scanB(2:NB,2)'-scanB(1:(NB-1),2)';
denom=sqrt(distX.^2+distY.^2); % the distances between every consecutive points [1*(N-1)]
Ax=scanA(:,1); % N*1
Ay=scanA(:,2);
Bx=scanB(:,1)'; % 1*N
By=scanB(:,2)';
distance=(distX(ones(NA,1),:).*(By(ones(NA,1),1:end-1)-Ay(:,ones(1,NB-1))) - ...
    (Bx(ones(NA,1),1:end-1)-Ax(:,ones(1,NB-1))).*distY(ones(NA,1),:))...
    ./denom(ones(NA,1),:);
% the distance is N*(N-1)
[mindist,index] = min(abs(distance)');% mindist is 1*N, index is 1*N 
indexes = [[1:NA]' index']; % N*2;
ind = find(mindist<thres^2);
indexes = indexes(ind,:);
%Computes the intersection point 
intscaB(:,1)=(-distY(index(ind))./denom(index(ind)))...
    .*distance(sub2ind([NA,NB-1],ind,index(ind)))+scanA(ind,1)';
intscaB(:,2)=(distX(index(ind))./denom(index(ind)))...
    .*distance(sub2ind([NA, NB-1],ind,index(ind)))+scanA(ind,2)';