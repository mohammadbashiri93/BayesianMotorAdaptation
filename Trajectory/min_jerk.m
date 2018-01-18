% Offline jerk with normalisation
% Can be converted to online (Page 6, Shadmeyer tutorial)


clc;
clear all;


% NOTE: everthing is represented in rows -> [x, y]

% Position Information
startPos  = [0, 0];
r_i = startPos;
targetPos = [0, 8];
r_f = targetPos;

% simulation params
t = 0;
dt = 0.01;    % [s]  -> step size
t_f = 1; % [s]
N = t_f/dt; 

r_list = zeros(N-1,2);
r1_list =zeros(N-1,2);
r2_list =zeros(N-1,2);
tau_list = zeros(N-1,1);

r  = startPos; % initial position
v  = [+0.0, 0.0]; % initial speed

r1_init = [2.2 0.0];
r2_init = [20.0 0.0];


r_past = [0 0];
r1_past = [2.2 0.0];
r1 = [2.2 0.0];
r = [0 0];
r2 = [20.0 0.0];

t = dt;

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
    a2 = D.*r2_init./2;
    a3 = (-(3*D^2)/2).*r2_init - 6*D.*r1_init + 10.*(targetPos - startPos);
    a4 = ((3*D^2)/2).*r2_init + 8*D.*r1_init - 15.*(targetPos - startPos);
    a5 = (-(D^2)/2).*r2_init - 3*D.*r1_init + 6.*(targetPos - startPos);
     
%     T = t_f;
%     tau = t; %dt;  % we are applying it every step, so the difference is dt
%     
%     % perturbation addition 
%     T_mat = [T^3 T^4 T^5; 3*T^2 4*T^3 5*T^4; 6*T 12*T^2 20*T^3]; %3x3
%     B = [(r_f-r_i-r1_init*T-(r2_init/2)*T^2); (-r1_init-r2_init*T); -r2_init]; % 3x2   I have to feed current state values
%     A = T_mat\B; %3x2
%       
%     a0 = r_i;
%     a1 = r1_init;
%     a2 = r2_init/2;
%     a3 = A(1,:);
%     a4 = A(2,:);
%     a5 = A(3,:);
%     
    
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


figure;
subplot(1,3,1);
plot(r_list(:,1), r_list(:,2));
subplot(1,3,2);
plot(r1_list(:,1), r1_list(:,2));
subplot(1,3,3);
plot(r2_list(:,1), r2_list(:,2));

figure;
x = vecnorm(r1_list,2,2);  % velocity norm
plot(x);
hold on;
y = r_list(:,2);% y component of position
plot(y);
hold on;
z = r_list(:,1); % x component of position
plot(z);
hold on;
zz = vecnorm(r2_list,2,2);
plot(zz);



