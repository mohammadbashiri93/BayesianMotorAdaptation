function posteriorSpace = UpdateBelief( priorSpace, currentState, Fapplied, Generalization )
%UPDATEBELIEF Summary of this function goes here
%   Detailed explanation goes here

%% Define global and other constant variables

global Vx; global Vy; global Fres;

posteriorSpace = priorSpace; 
likelihood = compLikelihood(Fapplied);

%% Update current state in with/without generalization

GF = compGF(currentState);

if (Generalization == true)
    strtVx = 1; 
    endVx = length(Vx);
    strtVy = 1;
    endVy = length(Vy);
else
    [indX, indY] = findStateInd(currentState);
    strtVx = indX;
    endVx = indX;
    strtVy = indY;
    endVy = indY;
end

for indX = strtVx:endVx
    for indY = strtVy:endVy
  
        Gcoef = GF(indX, indY);  
        
        stateLikelihood = likelihood * Gcoef;

        priorBelief = buildPriorBelief(indX, indY, priorSpace);

        posteriorBelief = stateLikelihood .* priorBelief;           
        pState = sum(posteriorBelief(:)) * Fres^2;            
        posteriorBelief = posteriorBelief/pState;

        pFx = sum(posteriorBelief, 1);
        Fx_mu = mean(pFx);
        Fx_sigma = std(pFx);

        pFy = sum(posteriorBelief, 2);
        Fy_mu = mean(pFy);
        Fy_sigma = std(pFy);

        posteriorSpace(1,1,indX,indY) = Fx_mu;
        posteriorSpace(1,2,indX,indY) = Fx_sigma;
        posteriorSpace(2,1,indX,indY) = Fy_mu;
        posteriorSpace(2,2,indX,indY) = Fy_sigma;
    end
end

end

