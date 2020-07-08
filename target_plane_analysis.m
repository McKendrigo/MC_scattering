% Target analysis
%%%%%%%%%%%%%%%%%

% This code analyses the coordinates of photons that have intercepted at
% target plane (outputs from fast_plane_analysis.m) and:

% a) Calculates the sum of packet weights that "hit" a square detector with
% a user-specified area, that is parallel to the target plane.

% b) Calculates the distribution of packet weights across a user-specified 
% area within the target plane.

% 1) Setting up
%%%%%%%%%%%%%%%

prompt = {'Target centre coordinates:','Target width (m):','Grid width (m):','Grid step size (m):'};
dlgtitle = 'Inputs required';
dims = [1 40];
definput = {'[0,0]','1e-3','1','0.05'};
trgt_plane_inputs = inputdlg(prompt,dlgtitle,dims,definput);

centre = str2num(trgt_plane_inputs{1}); % Target centre coordinates
max_deviation = 0.5 * str2num(trgt_plane_inputs{2}); % Maximum x/y deviation from the centre

x_min = centre(:,1) - max_deviation;
x_max = centre(:,1) + max_deviation;
y_min = centre(:,2) - max_deviation;
y_max = centre(:,2) + max_deviation;

% 2) Find all packets within target area and sum their packet weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target_hits = find(coordinates(:,1) >= x_min & coordinates(:,1) <= x_max...
    & coordinates(:,2) >= y_min & coordinates(:,2) <= y_max);
Rx_received_total = sum(weights(target_hits)); % This is the total sum of packets hitting the receiver.

% 3) Calculated distribution of packets within the specified 'grid'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_width = str2double(trgt_plane_inputs{3}); % Width of grid (m)
gridSize = str2double(trgt_plane_inputs{4}); % Grid step size (m)

% Define target area
targetLoc = [-0.5*grid_width, -0.5*grid_width, grid_width, grid_width];
%[x (lower left corner), y (lower left corner), width, height], (m)

% Create target grid
xGrid = targetLoc(1) + gridSize*(0:floor(targetLoc(3)/gridSize));
yGrid = targetLoc(2) + gridSize*(0:floor(targetLoc(4)/gridSize));

% Find packets within grid boundaries
grid_hits = find(abs(coordinates(:,1)) <= 0.5*grid_width & ... 
    abs(coordinates(:,2)) <= 0.5*grid_width); 

% Compute number of photon "hits" within each grid box
[nHits,~,~,binX,binY] = histcounts2(coordinates(grid_hits,1), coordinates(grid_hits,2), xGrid, yGrid); 

% Sum weights within each bin
[~, hitGroups, hitGroupID] = unique([binX,binY],'rows','stable');
totWeights = splitapply(@sum,weights(grid_hits),hitGroupID); 
ind = sub2ind(size(nHits),binX(hitGroups), binY(hitGroups));
weightMatrix = nHits; 
weightMatrix(ind) = totWeights;

% Add weighted hit plot
ax2 = subplot(1,1,1); 
I = imagesc(ax2,xGrid(1:end-1)+gridSize/2, yGrid(1:end-1)+gridSize/2,weightMatrix');
ax2.YDir = 'normal';
linkprop([ax,ax2],{'xlim','ylim','xtick','ytick','XTickLabel','YTickLabel'})
grid(ax2,'on')
axis(ax2,'equal')
axis(ax2,'tight')
cb2 = colorbar(ax2);
ax2.CLim = [0,max(totWeights)]; 
ax2.Colormap(1,:) = [1,1,1]; % This sets 0-values to white
ylabel(cb2,'Weight sum')
title(ax2,'Weighted hits')
text(ax2, xGridMat(hitIdx)+gridSize/2, yGridMat(hitIdx)+gridSize/2, compose('%.2f',weightMatrix(hitIdx)), ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'middle','Fontsize', 10, 'Color', 'r')
text(ax2, min(xlim(ax2)), min(ylim(ax2)), 'Numbers show sum of weights', 'VerticalAlignment', 'bottom')

% 
% [xq,yq] = meshgrid(-grid_width/2:grid_step:grid_width/2,...
%     -grid_width/2:grid_step:grid_width/2); % Target mesh grid
% 
% vq = griddata(coordinates(:,1),coordinates(:,2),weights,xq,yq);
% % Interpolate the 3D scattered data (2D coordinate and weight)
% 
% % Plot received intensity profile at the target plane.
% contourf(xq,yq,vq);
% xlabel('X (m)'); ylabel('Y (m)');
% set(gca, 'FontSize', 20);
% title('Interpolated received intensity at target');
