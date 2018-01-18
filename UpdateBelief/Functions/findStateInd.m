function [ indVx, indVy ] = findStateInd( currentState )
% FINDSTATEIND finds the closest state variable to the given 
% state (velocity) values

% INPUT:
%   - State variable (current state values)
%   - State variable range (Vx, Vy)
%
% OUTPUT
%   - CurrentVx, currentVy

global Vx;
global Vy;

n = 1; % number of decimal points

valVx = currentState(1);
valVy = currentState(2);

indVx = find(round(Vx,n) == round(valVx,n));
indVy = find(round(Vy,n) == round(valVy,n));

end

