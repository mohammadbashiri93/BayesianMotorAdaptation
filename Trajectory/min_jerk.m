% Offline jerk with normalisation
% Can be converted to online (Page 6, Shadmeyer tutorial)


clc;
clear all;


% NOTE: everthing is represented in rows -> [x, y]

% Position Information
startPos  = [0, 0];
targetPos = [0, 8];


% simulation params
t = 0;
dt = 0.1;    % [s]  -> step size
t_f = 5; % [s]
N = t_f/dt; 

r_list = zeros(N-1,2);
r1_list =zeros(N-1,2);
r2_list =zeros(N-1,2);
tau_list = zeros(N-1,1);

r  = startPos; % initial position
v  = [0, 0]; % initial speed

r1_init = [0 0];
r2_init = [0 0];


r_past = [0 0];
r1_past = [0 0];
r1 = [0 0];
r = [0 0];
r2 = [0 0];



% here is were the loop should be implemented
figure; grid; hold on

for i = 1:N-1

    % Implementation of minimum jerk methof
    D =  t_f;
    tau = t/D;  % we are applying it every step, so the difference is dt
     
    % some extra parameters which might be of interest
    q = [r' r1' r2']; % 2x3
    A = [0 0 -60/D^3; 1 0 -36/D^2; 0 1 -9/D]; % 3x3
    B = [0 0 60/D^3]; % 1x3
    B = repmat(B,2,1); %2x3
    q1 = q * A + B.* repmat(targetPos',1,3);  % q1 =[r1', r2', r3']   2x3
    
    a0 = startPos;
    a1 = D.*r1_init;
    a2 = D.*r2_init/2;
    a3 = (-(3*D^2)/2)*r2_init - 6*D*r1_init + 10*(targetPos - startPos);
    a4 = ((3*D^2)/2)*r2_init + 8*D*r1_init - 15*(targetPos - startPos);
    a5 = (-(D^2)/2)*r2_init - 3*D*r1_init + 6*(targetPos - startPos);
    
    a = [a0' a1' a2' a3' a4' a5']; %2x6
    
    taus = [1; tau; tau^2; tau^3; tau^4; tau^5]; % 6x1
    
    % update position
    r = a*taus;
    r = r';
    
    %r1 = [q1(1,1) q1(2,1)]; % this is the one used in this step
    %r2 = [q1(1,2) q1(2,2)]; % this is the one used in this step
    
    r1 = (r - r_past)/dt;
    r2 = (r1 - r1_past)/dt;
       
    % draw the trajectory
    scatter(r(1), r(2), '.'); xlim([-20,20]); ylim([-20,20]);
    % draw the start and target points
    scatter(startPos(1), startPos(2), 'filled', 'k');
    scatter(targetPos(1), targetPos(2), 'filled', 'r');     

    r_past = r;
    r1_past = r1;
    r_list(i,:) = r;
    r1_list(i,:) = r1;
    r2_list(i,:) = r2;
    tau_list(i) = tau;
    t =  t +dt;
    
    pause(0.001)
end