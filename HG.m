% DUV scattering function %
% See Ding et al., IEEE Journal on Selected Areas in Comms., 2009 %
function [PDFfit] = HG(ks_ray,ks_mie,upsilon,g,f)

theta = 0:(pi/1000):pi; % Scatting angles (Radians)
THETA = transpose(theta);

% Fit parameters, taken from Ding et al. %
mu = transpose(cos(theta));

% Rayleigh and Mie scattering coefficients, taken from Ding et al. %

ks = ks_ray + ks_mie;

% Set up the equations that describe Rayleigh and Mie scattering %
p_ray_numerator = 3*(1+(3*upsilon)+((1-upsilon)*mu.^2));
p_ray_denominator = 16*pi*(1+2*upsilon);

p_ray = p_ray_numerator/p_ray_denominator;

p_mie_factor = (1-g^2) / (4*pi);
p_mie_denom1 = (1+(g^2)-(2*g*mu)).^(3/2);
p_mie_denom2 = (1+g^2)^(3/2);
p_mie_num2 = 0.5*((3*mu.^2)-1);

p_mie = p_mie_factor*((1./p_mie_denom1) + f*((p_mie_num2./p_mie_denom2)));

p = ((ks_ray/ks)*p_ray) + ((ks_mie/ks)*p_mie); % This the overall scattering phase function

%%% Fit a linear spline to p vs. mu %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fitobj = fit(THETA,p.*sin(THETA),'linearinterp');
T = 2*pi*(integrate(fitobj,theta,0));

PDFfit = fit(transpose(T),transpose(theta),'linearinterp'); % Fit the Theta vs. T data
end