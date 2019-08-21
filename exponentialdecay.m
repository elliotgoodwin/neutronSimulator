function [ s, MFP_gradient, n, x, f, MFP_error] = exponentialdecay( a, MFP)
% A function that generates random numbers and samples them according to an
% exponential decay distribution

% first generate a table of uniformly distributed random numbers
u = rand(a, 1);

% mulitply uniformly distributed variable u by inverse CDF to get a
% distribution according to the original PDF
% MFP = mean free path 
% s = sampled value
s = -MFP.*log(u);

end

