%% We'll start by simulating a point mass in 2D and apply a force on it and see the trajectory
%%%%%%%%%%%%%%%%%%%%%% Minimum Jerk Method %%%%%%%%%%%%%%%%%%%%%%%
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
t_f = 50; % [s]
N      = t_f/dt; 

r  = startPos; % initial position
v  = [0, 0]; % initial speed

r_past = 0;
r1_past = 0;

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

    % I think we need to map the force directly to velocity
    v = f *10;
    
    % Implementation of minimum jerk methof
    D =  t_f - t;
    tau = t/D;
    
    r1 = (r - r_past)/dt;
    r2 = (r1 - r1_past)/dt;
    
    q = [r r1 r2]; % 2x3
    A = [0 0 -60/D^3; 1 0 -36/D^2; 0 1 -9/D]; % 3x3
    B = [0 0 60/D^3]; % 1x3
    B = repmat(B,2,1); %2x3
    
    a0 = r;
    a1 = D.*r1;
    a2 = D.*r2/2;
    a3 = (-(3*D^2)/2)*r - 6*D*r1 + 10*(target_pos - r);
    a4 = ((3*D^2)/2)*r + 8*D*r1 - 15*(target_pos - r);
    a5 = (-(D^2)/2)*r - 3*D*r1 + 6*(target_pos - r);
    
    
    q1 = q * A + B.* repmat(r',1,3);  % repmat has to be used. See the paper

        
    
    % update position
    r = r + v.*dt; % + 0.5.*a.*dt^2;
    
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

    r_past = r;
    r1_past = r1;
    
    pause(0.001)
end