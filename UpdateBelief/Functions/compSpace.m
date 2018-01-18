function [ complete_sspace ] = compSpace( sspace )
% COMPSPACE computes the complete 4D state space (including the force distributions)

% INPUT:
%   - state space
%   - Fx vector
%   - Fy vector

% OUTPUT:
%   - complete state space, includnig the force distribution (NOT the mu and sigma)

global Fx;
global Fy;

Fsize = length(Fx);

complete_sspace = zeros(Fsize, Fsize, size(sspace,3), size(sspace,4));

for indVx = 1:size(sspace,3)
    for indVy = 1:size(sspace,4)
        
        F = Gauss2d( Fx, sspace(1, 1, indVx, indVy), sspace(1, 2, indVx, indVy),...
                     Fy, sspace(2, 1, indVx, indVy), sspace(2, 2, indVx, indVy));
        
        complete_sspace(:,:,indVx, indVy) = F;
        
    end
end
        
end

