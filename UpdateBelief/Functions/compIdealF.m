function [ ideal_F ] = compIdealF( func )
% COMPIDEALF computes the ideal force compensation values

% This function computes the ideal force compensation values, given the
% state (velocity)

% INPUT:
%   - rotation matrix function (this is dynamics-specific)
%   - Vx vector
%   - Vx vector

% OUTPUT: 
%   - ideal force distribution. dimension: (2, Vx, Vy);

global Vx;
global Vy;

ideal_F = zeros(2, length(Vx), length(Vy));

for indVx = 1:length(Vx)
    for indVy = 1: length(Vy)
   
        ideal_F(:, indVx, indVy) = func([Vx(indVx) Vy(indVy)]);
        
    end
end

end

