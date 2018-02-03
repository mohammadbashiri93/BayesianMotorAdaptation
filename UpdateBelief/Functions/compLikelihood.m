function likelihood = compLikelihood( Fapplied )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global Fx; global Fy; global Likely_sigmas;

FappliedX = Fapplied(1);
FappliedY = Fapplied(2);

likely_mus = [FappliedX, FappliedY];

% Let's say dim(Fx) = p and dim(Fy) = q
[FX, FY] = meshgrid(Fx, Fy);  % q x p. 
likelihood = mvnpdf([FX(:), FY(:)], likely_mus, Likely_sigmas);

%likelihood = reshape(likelihood, length(Fx), length(Fy))'; Incorrect
%(Pranshul)
%1st index should be Fy, 2nd should be Fx
likelihood = reshape(likelihood, length(Fy), length(Fx)); % q x p

end

