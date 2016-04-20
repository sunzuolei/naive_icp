%function [pos,indexes,err]=icp(rawscanA, rawscanB,  relpos, thres_trans, nit, interp, plotflag)
%Implements the ICP algorithm with no interpolation
%
%INPUTS:
% rawscanA, rawscanB vectors with laser scans
% relpos intial relative position
% thres_trans threshold on maximum translation
% nit number of iterations
% interp flag to perform interpolation
%
%OUTPUTS:
% pos relative position
%
% Implemented by Fabio Tozeto Ramos

function pos=icp(rawscanA, rawscanB,  relpos_ba, thres_trans, nit)
pointsA = get_laser_points(rawscanA, 5.8); % [x; y]; (2*361)
pointsB = get_laser_points(rawscanB, 5.8);
pos = relpos_ba;


prevpos=relpos_ba;
for i=1:nit

    pointsA_rel_B=compound(pos,pointsA);
    indexes=nn_mutual(pointsA_rel_B',pointsB',thres_trans);


    
    

    
    [post(1),post(2),post(3)]=rel_pose(pointsA_rel_B(:,indexes(:,1))',pointsB(:,indexes(:,2))');
   
    pos=pos+post';
    posdist=sqrt((pos(1)-prevpos(1))^2+(pos(2)-prevpos(2))^2);
   
    if posdist<1e-8
        break;
    else
        prevpos=pos;
    end

end


