%%% Scatter packets - ISOTROPIC SCATTERING - update direction cosines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dir,psi] = scatter_packets(dir,prev_psi,active_index)

% First, we create the new scattering angles for each packet
t1 = rand(length(dir),1); % Dice roll for Elevation/Scattering (Psi) angles - all packets
t2 = rand(length(dir),1); % NEW CODE = Dice roll for Azimuthal scattering angle - all packets

t1(active_index==0) = NaN; % t1 for inactive packets set to "NaN"
t2(active_index==0) = NaN; % t2 for inactive packets set to "NaN"

scatter_psi = acos(1-(2*t1)); %%% psi formula for isotropic scattering
psi = scatter_psi + prev_psi; % New psi angle of travel, relative to z-axis
theta = 2*pi*t2; % NEW CODE = update azimuthal scattering angles

% See Leathers Section 4.2 %
% For photons that are very close to the z-axis (nearly vertical), it is
% preferable to update their direction cosines as follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vert_index = find(abs(dir(active_index,3)) >= 0.99999); % Indices of active packets that are close to vertical
other_index = find(abs(dir(active_index,3)) < 0.99999); % Indices of all other active packets that are not close to vertical

dir(vert_index,1) = sin(psi(vert_index,1)) .* cos(theta(vert_index,1));
dir(vert_index,2) = sin(psi(vert_index,1)) .* sin(theta(vert_index,1));
dir(vert_index,3) = cos(psi(vert_index,1));

% For all other packets we update the direction cosines as follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update those index cosines
x = dir(other_index,1); % Buffer the old x dir values, so that when x dir is updated it doesn't mess up the calculation of new y dir values!

dir(other_index,1) = (((dir(other_index,1).*dir(other_index,3))./sqrt((1-dir(other_index,3).^2))) .* (sqrt(1-(cos(psi(other_index,1))).^2)).*cos(theta(other_index,1))) + ((-dir(other_index,2)./sqrt((1-dir(other_index,3).^2))) .* (sqrt(1-(cos(psi(other_index,1))).^2)).*sin(theta(other_index,1))) + (dir(other_index,1) .* cos(psi(other_index,1)));
dir(other_index,2) = (((dir(other_index,2).*dir(other_index,3))./sqrt((1-dir(other_index,3).^2))) .* (sqrt(1-(cos(psi(other_index,1))).^2)).*cos(theta(other_index,1))) + ((x./sqrt((1-dir(other_index,3).^2))) .* (sqrt(1-(cos(psi(other_index,1))).^2)).*sin(theta(other_index,1))) + (dir(other_index,2) .* cos(psi(other_index,1)));
dir(other_index,3) = ((-sqrt((1-dir(other_index,3).^2))) .* (sqrt(1-(cos(psi(other_index,1))).^2)).*cos(theta(other_index,1))) + 0 + (dir(other_index,3) .* cos(psi(other_index,1)));

end