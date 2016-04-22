%Given two scans, plot the resulting association given by ICP
%ploticpassoc(scanA, scanB)
function ploticpassoc(scanA, scanB)

%Run ICP twice, with two different thresholds

[pos,indexes,err]=icp(scanA, scanB,  [0;0;0], 5, 50, 0, 0);
[pos,indexes,err]=icp(scanA, scanB,  pos, 1, 50, 0, 0);

if size(scanA,1)~=2
    scanA = get_laser_points(scanA, 70);
    scanB = get_laser_points(scanB, 70);
else
    scanA = scanA;
    scanB = scanB;
end

%Plot associations
plot_assoc_gen(scanA,scanB,indexes(:,1),indexes(:,2));