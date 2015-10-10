%[w,Tx,Ty] = rel_pose_w(scanA,scanB,w)
%Compute the relative position of two scans based on the weighted least-square
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
% Implemented by Fabio Tozeto Ramos 14/11/06
function [Tx,Ty,w] = rel_pose(scanA,scanB,w)
N = length(scanA);
%x_dashA = mean(scanA(:,1));
%x_dashB = mean(scanB(:,1));
%y_dashA = mean(scanA(:,2));
%y_dashB = mean(scanB(:,2));
%SxAxB = sum((scanA(:,1)-x_dashA*ones(N,1)).*(scanB(:,1)-x_dashB*ones(N,1)),1);
%SyAyB = sum((scanA(:,2)-y_dashA*ones(N,1)).*(scanB(:,2)-y_dashB*ones(N,1)),1);
%SxAyB = sum((scanA(:,1)-x_dashA*ones(N,1)).*(scanB(:,2)-y_dashB*ones(N,1)),1);
%SyAxB = sum((scanA(:,2)-y_dashA*ones(N,1)).*(scanB(:,1)-x_dashB*ones(N,1)),1);
xA = sum(scanA(:,1).*w);
xB = sum(scanB(:,1).*w);
yA = sum(scanA(:,2).*w);
yB = sum(scanB(:,2).*w);
Sw = sum(w);

SxAxB = sum(scanA(:,1).*scanB(:,1).*w) - xA*xB/Sw;
SyAyB = sum(scanA(:,2).*scanB(:,2).*w) - yA*yB/Sw;
SxAyB = sum(scanA(:,1).*scanB(:,2).*w) - xA*yB/Sw;
SyAxB = sum(scanA(:,2).*scanB(:,1).*w) - yA*xB/Sw;

x_dashA = xA/Sw;
x_dashB = xB/Sw;
y_dashA = yA/Sw;
y_dashB = yB/Sw;


w = atan2((SxAyB-SyAxB),(SxAxB+SyAyB)); 
Tx = x_dashB - (x_dashA*cos(w) - y_dashA*sin(w));
Ty = y_dashB - (x_dashA*sin(w) + y_dashA*cos(w));
end

