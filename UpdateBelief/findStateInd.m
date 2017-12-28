function [ indVx, indVy ] = findStateInd( currentVx, currentVy, Vx, Vy )
% FINDSTATEIND finds the closest state variable to the given 
% state (velocity) values

% INPUT:
%   - State variable (current state values)
%   - State variable range (Vx, Vy)
%
% OUTPUT
%   - CurrentVx, currentVy

valVx = round(currentVx,1);
valVy = round(currentVy,1);

indVx = find(Vx == valVx);
indVy = find(Vy == valVy);

end

