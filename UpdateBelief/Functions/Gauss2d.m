function [ F ] = Gauss2d( x1, x1_mu, x1_sigma, x2, x2_mu, x2_sigma)
%GAUSS2D Summary of this function goes here
%   Detailed explanation goes here

mu = [x1_mu x2_mu];

% covaraince matrix
Sigma = [x1_sigma 0; 0 x2_sigma]; % variance within the data, and variance with other

% velocity range (min:step:max)
% x1 = -3:.2:3; x2 = -3:.2:3;
[X1, X2] = meshgrid(x1, x2);

% computer the 
F = mvnpdf([X1(:) X2(:)], mu, Sigma);

% plot it
F = reshape(F,length(x2),length(x1)); 

% normalize (make it a probability distribution)
F = F ./ sum(F(:));

end

