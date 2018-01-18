function likelihood = compLikelihood( Fapplied )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global Fx; global Fy; global likely_sigmas;

FappliedX = Fapplied(1);
FappliedY = Fapplied(2);

likely_mus = [FappliedX, FappliedY];

[FX, FY] = meshgrid(Fx, Fy);
likelihood = mvnpdf([FX(:), FY(:)], likely_mus, likely_sigmas);

likelihood = reshape(likelihood, length(Fy), length(Fx));

end

