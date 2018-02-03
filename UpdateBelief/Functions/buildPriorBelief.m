function [ priorBelief ] = buildPriorBelief( indx, indy, priorSpace )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

global Fx; global Fy;

Fx_mus    = squeeze(priorSpace(1, 1, indx, indy));
Fx_sigmas = squeeze(priorSpace(1, 2, indx, indy));
Fy_mus    = squeeze(priorSpace(2, 1, indx, indy));
Fy_sigmas = squeeze(priorSpace(2, 2, indx, indy));

F_mus = [Fx_mus, Fy_mus];
F_sigmas = [Fx_sigmas 0;0 Fy_sigmas];

% Let's say dim(Fx) = p and dim(Fy) = q
[FX, FY] = meshgrid(Fx, Fy); % q x p

priorBelief = mvnpdf([FX(:), FY(:)], F_mus, F_sigmas);
%priorBelief = reshape(priorBelief, length(Fy), length(Fx))'; Incorrect
%(pramhul)

priorBelief = reshape(priorBelief, length(Fy), length(Fx)); % q x p
end

