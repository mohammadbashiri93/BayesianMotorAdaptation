function [ Fadapt ] = useBelief( currentState, sspace )
% USEBELIEF returns the adaptation force (the force applied by the subject)
%
% INPUTS: 
%   - Current state (current Vx, current Vy)
%   - Range of states (Vx, Vy), State space
%   - State space (sspace)

% OUTPUT: 
%   - Adaptation Force

% steps:
%
% 1) Find the force distribution for the given state (Vx, Vy)
% 2) Sample from this distribution to find Fx and Fy
% 3) 

% TODO: There is a need to implement some more code to deal with edge
% cases here
[indVx, indVy] = findStateInd(currentState);

% The mu and sigma for the force distribution of current state
mu    = [sspace(1,1,indVx, indVy) sspace(2,1,indVx, indVy)];
sigma = [sspace(1,2,indVx, indVy) 0;0 sspace(2,2,indVx, indVy)];

Fadapt = mvnrnd(mu,sigma);

end

