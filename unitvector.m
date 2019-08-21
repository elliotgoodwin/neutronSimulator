function [ r ] = unitvector(a)
%Unit vector
%   a function that generates isotropic unit vectors, that are uniformly
% distributed over the surface of a sphere, by normalising a 3D normal
% distribution

% The 3D normal distribution is proportional to exp(- (x^2+y^2+z^2)). 
% Convert to polar coordinates ==> exp(-(r^2)). The density is 
% only a function of the radius, not the angle, which means the points are 
% uniformly distributed among all angles

% generate normally distributed random numbers
r = randn(a,3);
% normalise
r = bsxfun(@rdivide, r, sqrt(sum(r.^2,2)));

end

