% Compute the Hessian of the error function for the IDC algorithm
function H = hessian_idc(scanA, scanB, Z, dx, dy, dtheta)
H = zeros(3);

%d2E/dx2
F1 = calc_E_linear(scanA,scanB,[Z(1)+2*dx Z(2) Z(3)]);
F2 = calc_E_linear(scanA,scanB,[Z(1) Z(2) Z(3)]);
F3 = F2;
F4 = calc_E_linear(scanA,scanB,[Z(1)-2*dx Z(2) Z(3)]);
H(1,1) = derive2(F1,F2,F3,F4,dx,dx);

%d2E/dy2
F1 = calc_E_linear(scanA,scanB,[Z(1) Z(2)+2*dy Z(3)]);
F4 = calc_E_linear(scanA,scanB,[Z(1) Z(2)-2*dy Z(3)]);
H(2,2) = derive2(F1,F2,F3,F4,dy,dy);

%d2E/dtheta2
F1 = calc_E_linear(scanA,scanB,[Z(1) Z(2) Z(3)+2*dtheta]);
F4 = calc_E_linear(scanA,scanB,[Z(1) Z(2) Z(3)-2*dtheta]);
H(3,3) = derive2(F1,F2,F3,F4,dtheta,dtheta);

%d2E/dxdy
F1 = calc_E_linear(scanA,scanB,[Z(1)+dx Z(2)+dy Z(3)]);
F2 = calc_E_linear(scanA,scanB,[Z(1)-dx Z(2)+dy Z(3)]);
F3 = calc_E_linear(scanA,scanB,[Z(1)+dx Z(2)-dy Z(3)]);
F4 = calc_E_linear(scanA,scanB,[Z(1)-dx Z(2)-dy Z(3)]);
H(1,2) = derive2(F1,F2,F3,F4,dx,dy);
H(2,1) = H(1,2);

%d2E/dthetadx
F1 = calc_E_linear(scanA,scanB,[Z(1)+dx Z(2) Z(3)+dtheta]);
F2 = calc_E_linear(scanA,scanB,[Z(1)+dx Z(2) Z(3)-dtheta]);
F3 = calc_E_linear(scanA,scanB,[Z(1)-dx Z(2) Z(3)+dtheta]);
F4 = calc_E_linear(scanA,scanB,[Z(1)-dx Z(2) Z(3)-dtheta]);
H(3,1) = derive2(F1,F2,F3,F4,dx,dtheta);
H(1,3) = H(3,1);

%d2E/dthetady
F1 = calc_E_linear(scanA,scanB,[Z(1) Z(2)+dy Z(3)+dtheta]);
F2 = calc_E_linear(scanA,scanB,[Z(1) Z(2)+dy Z(3)-dtheta]);
F3 = calc_E_linear(scanA,scanB,[Z(1) Z(2)-dy Z(3)+dtheta]);
F4 = calc_E_linear(scanA,scanB,[Z(1) Z(2)-dy Z(3)-dtheta]);
H(3,2) = derive2(F1,F2,F3,F4,dy,dtheta);
H(2,3) = H(3,2);


function der = derive2(F1,F2,F3,F4,d12,d34)
der = (F1-F2)/(4*d12*d34)-(F3-F4)/(4*d12*d34);
