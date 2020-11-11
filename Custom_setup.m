% Define function for initial photon packet setup,based on an imported 
% custom intensity profile.
% See Khalighi et al. IEEE Photonics Journal 2017        

% Imported data should be formated as .csv, column 1 = angle (degress),
% column 2 = relative intensity (0 to 1).
% Final entries should always be 90degrees = 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [source_qfit,q] = Custom_setup(theta_degrees,profile)
theta = deg2rad(theta_degrees); % Convert degrees to radians
int_fit = fit(theta,profile,'linearinterp'); % Interpolate imported data
fun = @(x)(int_fit(x)/pi).*sin(x); % Fit function (intensity/solid angle)

q = zeros(length(theta),1); % Pre-allocate array for speed
index = 1;
while index < length(theta)+1
    q(index,:) = 2*pi*integral(fun,0,theta(index),'ArrayValued',true); % Integrate to find relationship between q and emission angle, theta
    index = index + 1;
end

q(end) = 1; % Constrain output so that q=1 at 90 degrees

source_qfit = fit(q,theta,'linearinterp'); % Fit to pass to the rest of the script