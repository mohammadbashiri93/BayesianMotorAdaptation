function [ sspace_image ] = showSpace( ideal_force, sspace )

% The correct way to do this is to have the ideal force distribution for
% each state. Then compute the probability of this ideal force
% compensation in the force distribution given the state! or the other way 
% around 

% we will implement this in th following section

% INPUTS: 
% - Current state sapce
% - Ideal force compensation for each state in the state space: 
%   dimension: (2, Vx, Vy)

% OUTPUT:
% - state space progression


Vsize = size(sspace, 3);
sspace_image = zeros(Vsize, Vsize);

for indVx = 1:Vsize
    for indVy = 1:Vsize
        
        mu    = [sspace(1,1,indVx, indVy), sspace(2,1,indVx, indVy)];
        Sigma = [sspace(1,2,indVx, indVy) 0; 0 sspace(2,2,indVx, indVy)];
        
%         disp('shpwspace')
%         disp(Sigma)
%         
        sspace_image(indVx, indVy) = mvnpdf(ideal_force(:,indVx, indVy)', mu, Sigma);
        
    end
end
end
