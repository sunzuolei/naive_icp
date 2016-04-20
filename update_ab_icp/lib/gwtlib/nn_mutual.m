%function [indexes] = nn_mutual(scanA, scanB, thres)
%Compute the nearst neighbour of each point in scan A from scan B and
%vice-versa.
%
%INPUTS:
% scanA a N x 2 matrix with scan values for A 
% scanB a N x 2 matrix with scan values for B
% thres a threshold on maximum distance
%
%OUTPUTS:
% a matrix Nmatches x 2 with the corresponce indexes of the two scans
%     indexes(:,1) - the index of point in ScanA which can find its correspondence in ScanB
%     indexes(:,2) - the index of point in ScanB which correspond the points indicated by indexes(:,1)
%
% Implemented by Fabio Tozeto Ramos 18/08/05

function [indexes] = nn_mutual(scanA, scanB, thres)
warning off all;
%if nargin<3, thres = 0.5; end

NA = size(scanA, 1);% NA size [N, 2]; N=361;
NB = size(scanB, 1);% NB size [N, 2];
XA2 = sum(scanA.^2,2);% XA2 size [N, 1]; x^2+y^2;
XB2 = sum(scanB.^2,2);% XB2 size [N, 1];
XB2T = XB2';% 1 by N;
distance = XA2(:,ones(1,NB))+XB2T(ones(1,NA),:)-2*scanA(:,1)*scanB(:,1)'-2*scanA(:,2)*scanB(:,2)';

[mindist,indexes] = min(distance,[],2);
% mindist is N*1, the minimum distance between the nth point in A and all of the points in B.
% indexes is N*1, the index of the points in B which is close to the nth point in A

j=1;
for i=1:length(indexes) % i is the index of point in A
    if length(find(indexes==indexes(i)))==1 % if only the point in B which is close to nth point in A, is only close to nth point in A
        indexes2(j,:)=[i, indexes(i)]; % [idx of point in A, idx of point in B which is close to this point in A]
        j=j+1;
    end
end

%indexes = [[1:NA]' indexes'];
ind = (mindist(indexes2(:,1))<thres^2); % check whether the minimum distance is less than the gate
% the result is false or true. NOTE: not 0 or 1!!!
indexes = indexes2(ind,:);

% the following codes are for testing the above line.
% c=[1,3;4,6;5,7];
% ind=[false false true];  c(ind,:) % that is OK.
% but if ind=[0 0 1]; % that is wrong! since the index should be larger than 0.
% NOTE: false ~= 0!!!
