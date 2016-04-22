% Calculate the Error function given a set of pair of points and the
% relative position
function [E, Nmatches] = calc_E(PA, PB, relpos)
PAglobal = transform_to_global(PA, relpos);
ind = nn(PAglobal', PB',2);
Nmatches = size(ind,1);
matchedA = PAglobal(:,ind(:,1));
matchedB = PB(:,ind(:,2));
XA2 = sum(matchedA.^2,1);
XB2 = sum(matchedB.^2,1);
E = sum(XA2 + XB2 - 2*matchedA(1,:).*matchedB(1,:) - 2*matchedA(2,:).*matchedB(2,:),2);