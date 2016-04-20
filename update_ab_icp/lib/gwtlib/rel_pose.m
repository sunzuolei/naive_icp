%function [Tx,Ty,w] = rel_pose(scanA,scanB)
%Compute the relative position of two scans based on the least-square
%solution.
%
%INPUTS:
% scanA = N x 2 matrix where N is the number of points
% scanB = N x 2 matrix where N is the number of points
% scanA and scanB need to have the same size
%
%OUTPUTS:
% w relative rotation
% Tx relative translation in X
% Ty relative translation in y
%
% Implemented by Fabio Tozeto Ramos 17/08/05
function [Tx,Ty,w] = rel_pose(scanA,scanB)
N = length(scanA);
M = length(scanB);

%x_dashA = mean(scanA(:,1));
%x_dashB = mean(scanB(:,1));
%y_dashA = mean(scanA(:,2));
%y_dashB = mean(scanB(:,2));
%SxAxB = sum((scanA(:,1)-x_dashA*ones(N,1)).*(scanB(:,1)-x_dashB*ones(N,1)),1);
%SyAyB = sum((scanA(:,2)-y_dashA*ones(N,1)).*(scanB(:,2)-y_dashB*ones(N,1)),1);
%SxAyB = sum((scanA(:,1)-x_dashA*ones(N,1)).*(scanB(:,2)-y_dashB*ones(N,1)),1);
%SyAxB = sum((scanA(:,2)-y_dashA*ones(N,1)).*(scanB(:,1)-x_dashB*ones(N,1)),1);
xA = sum(scanA(:,1)); % sum of the X coords in scan A
xB = sum(scanB(:,1));
yA = sum(scanA(:,2));
yB = sum(scanB(:,2));

SxAxB = sum(scanA(:,1).*scanB(:,1)) - xA*xB/N; % sum(XAn*XBm)- sum(XA)*sum(XB)/N
SyAyB = sum(scanA(:,2).*scanB(:,2)) - yA*yB/N;
SxAyB = sum(scanA(:,1).*scanB(:,2)) - xA*yB/N;
SyAxB = sum(scanA(:,2).*scanB(:,1)) - yA*xB/N;

x_dashA = xA/N;
x_dashB = xB/N;
y_dashA = yA/N;
y_dashB = yB/N;


w = atan2((SxAyB-SyAxB),(SxAxB+SyAyB)); 
Tx = x_dashB - (x_dashA*cos(w) - y_dashA*sin(w));
Ty = y_dashB - (x_dashA*sin(w) + y_dashA*cos(w));
end

