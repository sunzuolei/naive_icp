% Compute the covariance matrix
function [covar, covar2] = comp_cov(scanA, scanB, relpos)

[Emin, Nmat ] = calc_E_linear(scanA, scanB, relpos);
%s2 = eye(3); 
s2 = Emin/(Nmat-3-1)
dx = 0.025; dy = 0.025;
dtheta = deg2rad(1);
H = hessian_idc(scanA, scanB, relpos, dx, dy, dtheta)
covar = inv(0.5*H)*s2; 

