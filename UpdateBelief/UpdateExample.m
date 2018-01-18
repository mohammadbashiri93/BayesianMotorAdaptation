%%
% This code contains an example of updating the brior for every state
% we need to creat a 4D matrix which coontains every state
%
% This piece of code must be able to update and viualize the state space,
% as well as the force distribution of each state
%
% We are assuming Gaussian distributions, we will represent any
% Gaussian distribution using a mean and variance.So, our state space for 
% mean, as well as variance, would have 4 dimensions: (2, 2, Vx, Vy). 
% The "2" is because we have mean, and variance, for Fx and Fy (multivariate Gaussian)

%% Initializations

% initialize the global variables (mainly contains the biological aspects 
% of the subject)
InitGlobals()

% We have some implicit memory saved! and we can retrieve that memory.
% it can be based on textual cue!
memory     = retrieve_memory(1);

% some experiment-specific parameters
experiment = Exp_params(1);


global Vx;
global Vy;
Vsize = length(Vx);

sspace = zeros(2, 2, Vsize, Vsize);
sspace(:,1,:,:) = memory.Fmus;
sspace(:,2,:,:) = memory.Fsigmas;
% let's change some of the sigma:
% sspace(:, 2, 30:70, 30:70) = 0.5;
sspace(:, 2, :, :) = 0.2;

% Compute ideal force field values for each state
idealF = compIdealF(experiment.compF); 

%% Display the state space progression

sspace_image = showSpace(idealF, sspace);
figure; imagesc(Vx, Vy, sspace_image); xlabel('V_x'); ylabel('V_y'); axis xy
colorbar;

%% test useBelief function

currentState = [1.2, 1.2];
Fadapt = useBelief(currentState, sspace);

%% test updateBelief function


%% test the state space progression

sspace(:,1,5:10,5:10) = 0.9*idealF(:,5:10,5:10);  % 90% adapted
sspace_image = showSpace(idealF, sspace);
figure; imagesc(Vx, Vy, sspace_image); xlabel('V_x'); ylabel('V_y'); axis xy
colorbar;
