% Calculate the linearised Error function given a set of pair of points and the
% relative position
function [E, Nmatches] = calc_E_linear(PA, PB, relpos)
PAglobal = trans_global_linear(PA, relpos);
ind = nn(PAglobal', PB',1);
Nmatches = size(ind,1);
matchedA = PAglobal(:,ind(:,1));
matchedB = PB(:,ind(:,2));
XA2 = sum(matchedA.^2,1);
XB2 = sum(matchedB.^2,1);
E = sum(XA2 + XB2 - 2*matchedA(1,:).*matchedB(1,:) - 2*matchedA(2,:).*matchedB(2,:),2);

function p = trans_global_linear (p, relpos)
rot = [1 -relpos(3); relpos(3) 1];
p(1:2,:) = rot*p(1:2,:);

% translate
p(1,:) = p(1,:) + relpos(1);
p(2,:) = p(2,:) + relpos(2);

% if p is a pose and not a point
if size(p,1)==3
   p(3,:) = pi_to_pi(p(3,:) + relpos(3));
end
