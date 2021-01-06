% Target analysis
%%%%%%%%%%%%%%%%%

% This code analyses the coordinates of photons that have intercepted at
% target plane (outputs from fast_plane_analysis.m) and:

% a) Calculates the sum of packet weights that "hit" a square detector with
% a user-specified area, that is parallel to the target plane.

% b) Calculates the distribution of packet weights across a user-specified 
% area within the target plane.

% Info from:
% https://tinyurl.com/ycte3xbr

% 1) Setting up
%%%%%%%%%%%%%%%
function [Rx_received_total,xGrid,yGrid,weightMatrix] = target_plane_analysis(coordinates,hitweights,trgt_corners,trgt_plane,grid_num,Rx_previous,weightMatrix_previous,savetype,save_yn,outputfilename)

xmin = trgt_corners(1);
xmax = trgt_corners(2);
ymin = trgt_corners(3);
ymax = trgt_corners(4);
zmin = trgt_corners(5);
zmax = trgt_corners(6);

% Check in which plane the target lies, and find packets within the target
% area.
switch trgt_plane
    case 'x'
        % Define target area [(lower left corner),(lower left corner), width, height], (m)
        targetLoc = [ymin, zmin, ymax-ymin, zmax-zmin];
        % Find packets within grid boundaries
        grid_hits = find(coordinates(:,2) >= ymin && coordinates(:,2) <= ymax...
            && coordinates(:,3) >= zmin && coordinates(:,3) <= zmax); 
    case 'y'
        targetLoc = [xmin, zmin, xmax-xmin, zmax-zmin];
        % Find packets within grid boundaries
        grid_hits = find(coordinates(:,1) >= xmin && coordinates(:,1) <= xmax...
            && coordinates(:,3) >= zmin && coordinates(:,3) <= zmax); 
    case 'z'
        targetLoc = [xmin, ymin, xmax-xmin, ymax-ymin];
        % Find packets within grid boundaries
        grid_hits = find(coordinates(:,1) >= xmin && coordinates(:,1) <= xmax...
            && coordinates(:,2) >= ymin && coordinates(:,2) <= ymax); 
end

% 2) Calculate the distribution of packets within the specified 'grid'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create target grid
xGrid = linspace(targetLoc(1),targetLoc(1)+targetLoc(3),grid_num+1);
yGrid = linspace(targetLoc(2),targetLoc(1)+targetLoc(4),grid_num+1);

% Compute number of photon "hits" within each grid box
switch trgt_plane
    case 'x'
        [nHits,~,~,binX,binY] = histcounts2(coordinates(grid_hits,2), coordinates(grid_hits,3), xGrid, yGrid);
    case 'y'
        [nHits,~,~,binX,binY] = histcounts2(coordinates(grid_hits,1), coordinates(grid_hits,3), xGrid, yGrid);
    case 'z'
        [nHits,~,~,binX,binY] = histcounts2(coordinates(grid_hits,1), coordinates(grid_hits,2), xGrid, yGrid);
end

% Check if there are previous results to add
if isempty(weightMatrix_previous) == 1
    weightMatrix_previous = zeros(length(xGrid)-1,length(yGrid)-1); % If so, set all values to 0
end

% Sum weights within each bin
[~, hitGroups, hitGroupID] = unique([binX,binY],'rows','stable');
totWeights = splitapply(@sum,hitweights(grid_hits),hitGroupID); 
ind = sub2ind(size(nHits),binX(hitGroups), binY(hitGroups));
weightMatrix = nHits; 
weightMatrix(ind) = totWeights;
weightMatrix = weightMatrix + weightMatrix_previous; % Add results from this loop to the previous loop

% 3) Find all packets within target area and sum their packet weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch trgt_plane
    case 'x'
        target_hits = find(coordinates(:,2) >= ymin && coordinates(:,2) <= ymax...
            && coordinates(:,3) >= zmin && coordinates(:,3) <= zmax);
    case 'y'
        target_hits = find(coordinates(:,1) >= xmin & coordinates(:,1) <= xmax...
            & coordinates(:,3) >= zmin & coordinates(:,3) <= zmax);
    case 'z'
        target_hits = find(coordinates(:,1) >= xmin & coordinates(:,1) <= xmax...
            & coordinates(:,2) >= ymin & coordinates(:,2) <= ymax);
end

Rx_received_total = sum(hitweights(target_hits)); % This is the total sum of packets hitting the receiver
Rx_received_total = Rx_received_total + Rx_previous; % Add the previous total to the new total

% Save results (if this is the final loop to analyse)
if save_yn == "yes"
    if savetype == "Results only"
    save(outputfilename,'xGrid','yGrid', ...
    'Rx_received_total','weightMatrix'); % Save output
    elseif savetype == "All"
    save(outputfilename,'coordinates','hitweights','xGrid','yGrid', ...
    'Rx_received_total','weightMatrix'); % Save output
    end
end


