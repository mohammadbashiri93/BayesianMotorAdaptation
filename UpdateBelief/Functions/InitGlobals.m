function [ ] = InitGlobals()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global Vmin; Vmin = -100;
global Vmax; Vmax = 100;
global Vres; Vres = 1;

global Vx; Vx = Vmin:Vres:Vmax;
global Vy; Vy = Vx;

global Fmin; Fmin = -100;
global Fmax; Fmax = 100;
global Fres; Fres = 1;

global Fx; Fx = Fmin:Fres:Fmax;
global Fy; Fy = Fx;

sensory_sigma = 10^2;
gen_sigma = 5^2;
global Likely_sigmas; Likely_sigmas = [sensory_sigma 0; 0 sensory_sigma];
global GF_sigmas; GF_sigmas = [gen_sigma 0; 0 gen_sigma];

global m; m = 0.150; % hand mass [kg]

end

