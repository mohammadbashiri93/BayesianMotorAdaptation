% in this code we'll try to expand the single state simulation to
% trajectories.

% let's start by showing the single state adaptation


close all

parad = [zeros(1,50) ones(1,500)  -1.*ones(1,50) zeros(1,100) ones(1,120)];
% parad = [zeros(1,50) linspace(0,1,50)  -1.*ones(1,14) zeros(1,50) ones(1,120)];

V     = 0.999;
Alpha = 0.01;
x     = zeros(1,numel(parad));

% single state
for i = 2:numel(parad)
    error = parad(i-1) - x(i-1) ;
    x(i) = V*x(i-1) + Alpha*error;
end
figure; hold on;
plot(parad, 'k');
plot(x);
legend('Paradigm', 'Adaptation', 'Location', 'southwest')
% title('single state')

%% What we would like to have:

% 1) defined force field field
% 2) represent hand as a point mass
% And then having the hand adapt to the force field

%% We'll start by simulating a point mass in 2D and apply a force on it and see the trajectory

% NOTE: everthing is represented in rows -> [x, y]

% define mass ad force
m  = 0.250;          % [kg] -> average mass of the hand
f  = [0, 0.01];   % [N]  -> force applied on the hand


% force field
angle = 90;

% simulation params
dt     = 0.1;    % [s]  -> step size
period = 40; % [s]
N      = period/dt; 

r  = [0, 0]; % initial position
v  = [0, 0]; % initial speed

% here is were the loop should be implemented
figure; grid; hold on
for i = 1:N
    
    % calculate the acceleration
    a = f./m;
    
    % update velocity
    v = a.*dt + v;
    
    % update position
    r = r + v.*dt + 0.5.*a.*dt^2; % in principle, we should not have acceleration term here
    
    scatter(r(1), r(2)); xlim([-20,20]); ylim([-20,20]);
    
    pause(0.01)
    
    
    if mod(i,30)==0
        f = (rand(1,2)-0.5)./10; %(rand(1,2)-0.5)./10  % 0.07 .* rotateMat(v, angle);
    end
    
end

%% Now we will add a target

%% We'll start by simulating a point mass in 2D and apply a force on it and see the trajectory
% Two problems left: the task needs to be completed in the given time
% (modify force accordingly)
% make it adaptable to different force fields
%Solution: Implement the min jerk method. It automatically includes the
%restoring force. However, with perturbation, I would have to introduce
%deviation in position


% NOTE: everthing is represented in rows -> [x, y]

% define mass ad force
m  = 0.250;          % [kg] -> average mass of the hand
% f  = [0, 0.01];   % [N]  -> force applied on the hand


% force field
angle = 90;
type = 'NULL'; % 'FORCEFIELD'

% Position Information
startPos  = [0, 0];
targetPos = [0, 8];

% forces
f_goal        = [0 0];
f_hand        = [0 0];
f_forcefield  = [0 0];
f             = [0 0];

% simulation params
dt     = 0.1;    % [s]  -> step size
period = 50; % [s]
N      = period/dt; 

r  = startPos; % initial position
v  = [0, 0]; % initial speed

% here is were the loop should be implemented
figure; grid; hold on

for i = 1:N

    % compute total force applied on the hand    
    f_goal  = abs(r - startPos).* (targetPos - r) .* 0.01; % must be goal-dependent (HOW DO I CALCULATE THIS?)
    % if we keep it the above way, force would be max when most far from
    % the target. I does not happen like this in real experiments
    
    if (i==1)
        f_goal = (targetPos - r) .* 0.01;
    end
    
    f_adapt = -1 .* f_forcefield;
    f_hand  = f_adapt + f_goal; 
    f       = f_hand + f_forcefield;  % This is the net force influencing the motion
    
    % calculate the acceleration 
    %a = f./m;  % can't be done. It's not an external system under observation
    
    % update velocity
    %v = a.*dt + v;  % not correct
    
    % I think we need to map the force directly to velocity
    v = f *10;
    
    % update position
    r = r + v.*dt + 0.5.*a.*dt^2;
    
    % draw the trajectory
    scatter(r(1), r(2), '.'); xlim([-20,20]); ylim([-20,20]);
    % draw the start and target points
    scatter(startPos(1), startPos(2), 'filled', 'k');
    scatter(targetPos(1), targetPos(2), 'filled', 'r');
    % draw the force direction
    quiver(r(1), r(2), f_forcefield(1), f_forcefield(2),9);
      
    if mod(i,30)==0
        f_forcefield = 0.07 .* rotateMat(v, angle); %(rand(1,2)-0.5)./10
    end

    pause(0.001)
end