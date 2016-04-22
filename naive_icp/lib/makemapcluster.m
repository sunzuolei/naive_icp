%Uses ICP to create a map where every association of points is
%represented by a single point with position given by the mean
%of all instances.
%function [laserpts,assocs,poses,clpts] = makemapcluster(LSR)
function [laserpts,assocs,poses,clpts] = makemapcluster(LSR)

figure(1), axis equal, axis([-350 500 -450 450]); grid, hold on
laser=LSR.Range;
% Configurable values
INTERP = 0; % 0 - simple ICP, 1 - ICP with interpolation
NIT = 50;   % number of ICP iterations for each scan-match
GATE1 = 3; % 1st data association gate for each point in scan
GATE2 = 1;% 2nd data association gate for each point in scan
N = size(laser,1); % number of laser scans
laserpts = cell(1,N); % laser points
assocs = cell(1,N); % associations
poses = zeros(3,N); % poses
clpts = zeros(2,1e5); % clustered points
indpts = zeros(1,1e5); % indexes of the points in the current scan
mpts = zeros(1,1e5); % number of points in each cluster
ncl = 0; % index for the latest clustered point added

% Constants
MAXR = 70; % this laser data times-out at ~7.9 m

% Initialise
x = [0;0;0];
delta = [0;0;0];
p1=laser(1,:);
% Loop
for i=2:N
    p2 = laser(i,:); disp(i);
    %laserpts{i} = get_laser_points(p2, MAXR);
    delta = icp(p1, p2, delta, GATE1, 15,0,0);
    [delta, ind] = icp(p1, p2, delta, GATE2, NIT, INTERP,0);
    assocs{i} = ind;
    %delta = idc(p1, p2, delta, GATE2, NIT);
    
    x = transform_to_global(-delta, x);
    poses(:,i) = x;
    
    p1 = p2;

    % plots
    pt = get_laser_points(p2, MAXR);
    pg = transform_to_global(pt, x);
    laserpts{i} = pg;
    plot(x(1),x(2), '+', pg(1,:), pg(2,:), '.','erasemode','none');
    drawnow
end
figure(1); hold off

figure(2), axis equal, axis([-350 500 -450 450]); grid;

% Cluster the points
p1 = get_laser_points(laser(1,:), MAXR);


for i=2:N
    %p2 = get_laser_points(laser(i,:), MAXR);
    p2 = laserpts{i};
    cd1 = assocs{i}(:,1);% connected points
    cd2 = assocs{i}(:,2);
    %Fix indexes for points not associated
    indpts(~ismember(indpts,cd1))=0;
    
    % For each connection
    for j=1:length(cd1)    
        % Update the position if the point is already in the map
        if ismember(cd1(j),indpts) 
            ind = find(indpts==cd1(j));
            ind = ind(1); %hack if there is more than one point
            indpts(ind) = cd2(j);           
            clpts(1,ind) = (mpts(ind)*clpts(1,ind)+p2(1,cd2(j)))/(mpts(ind)+1);
            clpts(2,ind) = (mpts(ind)*clpts(2,ind)+p2(2,cd2(j)))/(mpts(ind)+1);
            mpts(ind) = mpts(ind)+1;
            
        %Add points if connected with previous scan    
        %elseif ismember(cd2(j),indpts)
        %    ind = find(indpts==cd2(j));
        %    ind = ind(1);
        %    clpts(1,ind) = (mpts(ind)*clpts(1,ind)+p1(1,cd1(j)))/(mpts(ind)+1);
        %    clpts(2,ind) = (mpts(ind)*clpts(2,ind)+p1(2,cd1(j)))/(mpts(ind)+1);
        %    mpts(ind) = mpts(ind)+1;
        else
            %Include a new point in the map
            ncl=ncl+1;
            mpts(ncl)=2;
            clpts(1,ncl) = (p1(1,cd1(j))+p2(1,cd2(j)))/2;
            clpts(2,ncl) = (p1(2,cd1(j))+p2(2,cd2(j)))/2;
            indpts(ncl) = cd2(j);
        end                
    end
    p1 = p2; 
end

plot(poses(1,:),poses(2,:),'b+',clpts(1,:),clpts(2,:),'g.');



