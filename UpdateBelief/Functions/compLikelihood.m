function likelihood = compLikelihood( Fapplied )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global Fx; global Fy; global Likely_sigmas;

FappliedX = Fapplied(1);
FappliedY = Fapplied(2);

likely_mus = [FappliedX, FappliedY];

[FX, FY] = meshgrid(Fx, Fy);
likelihood = mvnpdf([FX(:), FY(:)], likely_mus, Likely_sigmas);

likelihood = reshape(likelihood, length(Fx), length(Fy))';

end

