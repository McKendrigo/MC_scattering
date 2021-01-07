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
function [Rx_received_total,xGrid,yGrid,zGrid,weightMatrix] = target_plane_analysis(coordinates,hitweights,trgt_corners,trgt_plane,grid_num,Rx_previous,weightMatrix_previous,savetype,save_yn,outputfilename)

   
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
        grid_hits = find(coordinates(:,2) >= ymin & coordinates(:,2) <= ymax...
            & coordinates(:,3) >= zmin & coordinates(:,3) <= zmax);
        % Create target grid
        xGrid = [];
        yGrid = linspace(targetLoc(1),targetLoc(1)+targetLoc(3),grid_num+1);
        zGrid = linspace(targetLoc(2),targetLoc(2)+targetLoc(4),grid_num+1);
    case 'y'
        targetLoc = [xmin, zmin, xmax-xmin, zmax-zmin];
        % Find packets within grid boundaries
        grid_hits = find(coordinates(:,1) >= xmin & coordinates(:,1) <= xmax...
            & coordinates(:,3) >= zmin & coordinates(:,3) <= zmax);
        % Create target grid
        xGrid = linspace(targetLoc(1),targetLoc(1)+targetLoc(3),grid_num+1);
        yGrid = [];
        zGrid = linspace(targetLoc(2),targetLoc(2)+targetLoc(4),grid_num+1);
    case 'z'
        targetLoc = [xmin, ymin, xmax-xmin, ymax-ymin];
        % Find packets within grid boundaries
        grid_hits = find(coordinates(:,1) >= xmin & coordinates(:,1) <= xmax...
            & coordinates(:,2) >= ymin & coordinates(:,2) <= ymax);
        % Create target grid
        xGrid = linspace(targetLoc(1),targetLoc(1)+targetLoc(3),grid_num+1);
        yGrid = linspace(targetLoc(2),targetLoc(2)+targetLoc(4),grid_num+1);
        zGrid = [];
end

if isempty(grid_hits) == 0
    % 2) Calculate the distribution of packets within the specified 'grid'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    % Compute number of photon "hits" within each grid box
    switch trgt_plane
        case 'x'
            [nHits,~,~,binY,binZ] = histcounts2(coordinates(grid_hits,2), coordinates(grid_hits,3), yGrid, zGrid);
            % Sum weights within each bin
            [~, hitGroups, hitGroupID] = unique([binY,binZ],'rows','stable');
            totWeights = splitapply(@sum,hitweights(grid_hits),hitGroupID); 
            ind = sub2ind(size(nHits),binY(hitGroups), binZ(hitGroups));
        case 'y'
            [nHits,~,~,binX,binZ] = histcounts2(coordinates(grid_hits,1), coordinates(grid_hits,3), xGrid, zGrid);
            % Sum weights within each bin
            [~, hitGroups, hitGroupID] = unique([binX,binZ],'rows','stable');
            totWeights = splitapply(@sum,hitweights(grid_hits),hitGroupID); 
            ind = sub2ind(size(nHits),binX(hitGroups), binZ(hitGroups));
        case 'z'
            [nHits,~,~,binX,binY] = histcounts2(coordinates(grid_hits,1), coordinates(grid_hits,2), xGrid, yGrid);
            % Sum weights within each bin
            [~, hitGroups, hitGroupID] = unique([binX,binY],'rows','stable');
            totWeights = splitapply(@sum,hitweights(grid_hits),hitGroupID); 
            ind = sub2ind(size(nHits),binX(hitGroups), binY(hitGroups));
    end

    % Check if there are previous results to add
    if isempty(weightMatrix_previous) == 1 % If not, set all values to 0
        switch trgt_plane
            case 'x'
                weightMatrix_previous = zeros(length(yGrid)-1,length(zGrid)-1); 
            case 'y'
                weightMatrix_previous = zeros(length(xGrid)-1,length(zGrid)-1); 
            case 'z'
                weightMatrix_previous = zeros(length(xGrid)-1,length(yGrid)-1); 
        end
    end

    weightMatrix = nHits; 
    weightMatrix(ind) = totWeights; % Add weights in each "bin" to the results
    weightMatrix = weightMatrix + weightMatrix_previous; % Add results from this loop to the previous loop

    % 3) Find all packets within target area and sum their packet weights
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch trgt_plane
        case 'x'
        target_hits = find(coordinates(:,2) >= ymin & coordinates(:,2) <= ymax...
            & coordinates(:,3) >= zmin & coordinates(:,3) <= zmax);
        case 'y'
        target_hits = find(coordinates(:,1) >= xmin & coordinates(:,1) <= xmax...
            & coordinates(:,3) >= zmin & coordinates(:,3) <= zmax);
        case 'z'
        target_hits = find(coordinates(:,1) >= xmin & coordinates(:,1) <= xmax...
            & coordinates(:,2) >= ymin & coordinates(:,2) <= ymax);
    end

    Rx_received_total = sum(hitweights(target_hits)); % This is the total sum of packets hitting the receiver
    Rx_received_total = Rx_received_total + Rx_previous; % Add the previous total to the new total
    
else
    
end

% Save results (if this is the final loop to analyse)
if save_yn == "yes"
    if savetype == "Results only"
    save(outputfilename,'xGrid','yGrid','zGrid','trgt_corners', ...
    'trgt_plane','Rx_received_total','weightMatrix'); % Save output
    elseif savetype == "All"
    save(outputfilename,'coordinates','hitweights','xGrid','yGrid', ...
    'zGrid','Rx_received_total','weightMatrix',...
    'trgt_corners','trgt_plane'); % Save output
    end
end


