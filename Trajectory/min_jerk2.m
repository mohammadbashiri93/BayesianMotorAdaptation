% offline min jerk (Page 4, shadmeyer tutorial)

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


r  = startPos; % initial position
v  = [0, 0]; % initial speed

r1 = [0 0];
r2 = [0 0];


r_past = [0 0];
r1_past = [0 0];



% here is were the loop should be implemented
figure; grid; hold on

for i = 1:N-1    
    % some extra parameters which might be of interest
    

    A = [1; (10*(t/t_f)^3 -15*(t/t_f)^4 + 6*(t/t_f)^5)];
    X = [startPos' (targetPos-startPos)'];

    r = X*A;
    r = r';
    
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
    t =  t +dt;
    pause(0.001)
end

%%  
x = r1_list(:,2);
y = r2_list(:,2);

figure;
subplot(1,2,1);
plot(x);
subplot(1,2,2);
plot(y);