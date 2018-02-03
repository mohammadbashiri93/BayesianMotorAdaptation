function [posteriorSpace] = UpdateBelief( priorSpace, currentState, Fapplied, Generalization )
%UPDATEBELIEF Summary of this function goes here
%   Detailed explanation goes here

%% Define global and other constant variables

global Vx; global Vy; global Fres;
global Fx; global Fy;
global Likely_sigmas;

posteriorSpace = priorSpace; 
%%%likelihood = compLikelihood(Fapplied);  % 1st index: Fy; 2nd index: Fx

%% Update current state in with/without generalization
%%% Three percentage comments have been made by Pranshul

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

%for indX = strtVx:endVx   Incorrect (Pranshul)
for indX = strtVx:endVx
%    for indY = strtVy:endVy  Incorrect (pranshul)
    for indY = strtVy:endVy
        Gcoef = GF(indX, indY); % 1st index: Vx; 2nd index: Vy
        
%%%        stateLikelihood = likelihood * Gcoef;  % TODO: replace with Gcoef.

%%%        priorBelief = buildPriorBelief(indX, indY, priorSpace); % 1st index: Fy; 2nd index: Fx 
%  
%         disp('in update: before')
%         priorSpace(:,1,indX,indY)
        
%%%        posteriorBelief = stateLikelihood * 10 .* priorBelief; %  1st index: Fy; 2nd index: Fx

%%%        pState = sum(posteriorBelief(:)) * Fres^2;            
%%%        posteriorBelief = posteriorBelief/pState;

        
        
%         figure; 
%         subplot(311); surf(priorBelief); shading interp; view(2); axis tight
%         subplot(312); surf(stateLikelihood); shading interp; view(2); axis tight
%         subplot(313); surf(posteriorBelief); shading interp; view(2); axis tight
%         


%%%         pFx      = sum(posteriorBelief, 2)';
%%%         pFx      = pFx/sum(pFx);
%%%         Fx_mu    = sum(Fx.*pFx);
%%%         Fx_ex2   = sum(Fx.^2.*pFx);
%%%         Fx_sigma = Fx_ex2 - Fx_mu;
         
 
%%%         pFy      = sum(posteriorBelief, 1);
%%%         pFy      = pFy/sum(pFy);
%%%         Fy_mu    = sum(Fy.*pFy);
%%%         Fy_ex2   = sum(Fy.^2.*pFy);
%%%         Fy_sigma = Fy_ex2 - Fy_mu;


%         figure;
%         subplot(121); plot(Fx, pFx);
%         subplot(122); plot(Fy, pFy);

%       Sensory information    (Pranshul added)
        Fx_mu1 = Fapplied(1);
        Fy_mu1 = Fapplied(2);
        Fx_std1 = Likely_sigmas(1,1);
        Fy_std1 = Likely_sigmas(2,2);
    
        
%       Prior information (Pranshul added)
        Fx_mu2 = squeeze(priorSpace(1, 1, indX, indY));
        Fy_mu2 = squeeze(priorSpace(2, 1, indX, indY));
        Fx_std2 = squeeze(priorSpace(1, 2, indX, indY));
        Fy_std2 = squeeze(priorSpace(2, 2, indX, indY));

%       calculation of new parameters: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.583.3007&rep=rep1&type=pdf
        Fx_mu = Gcoef* (Fx_mu1 * Fx_std2^2 + Fx_mu2 * Fx_std1^2)/(Fx_std1^2 + Fx_std2^2);
        Fy_mu = Gcoef* (Fy_mu1 * Fy_std2^2 + Fy_mu2 * Fy_std1^2)/(Fy_std1^2 + Fy_std2^2);
        Fx_sigma = Gcoef * sqrt((Fx_std1^2 * Fx_std2^2)/(Fx_std1^2 + Fx_std2^2));
        Fy_sigma = Gcoef * sqrt((Fy_std1^2 * Fy_std2^2)/(Fy_std1^2 + Fy_std2^2));
        
        if Fx_sigma < 0.2
            Fx_sigma = .2;
        end
        if Fy_sigma < 0.2
            Fy_sigma = .2;
        end

        posteriorSpace(1,1,indX,indY) = Fx_mu;
        posteriorSpace(1,2,indX,indY) = Fx_sigma;
        posteriorSpace(2,1,indX,indY) = Fy_mu;
        posteriorSpace(2,2,indX,indY) = Fy_sigma;
%         disp('in update: after')
%         posteriorSpace(:, 1, indX, indY)

    end
end

end

