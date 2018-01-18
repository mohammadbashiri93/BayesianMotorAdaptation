%online min jerk without normalisation
%derived equations myself with the help of shadmeyer tutorial 
%I think we need to have online as well offline method. Then we need to
%weight them properly somehow. This is because particpants will not simply
%start using a new jerk method. They will probably use it in combination of
%what they are already using

clc;
clear all;


% NOTE: everthing is represented in rows -> [x, y]

% Position Information
r_i  = [0, 0];
r_f = [0, 8];


% simulation params
t =0;
dt = 0.01;    % [s]  -> step size
t_f = 1; % [s]
N = round((t_f)/dt); 

r_off_list = zeros(N,2);
r_on_list = zeros(N,2);  % online
r_list = zeros(N,2);
r_list(N,:) = r_f;
r1_list =zeros(N,2);
r2_list =zeros(N,2);
tau_list = zeros(N,1);

r_on_lists = 8*ones(N,2,N+1);

r  = r_i; % initial position
v  = [0, 0]; % initial speed

r1_fix = [0. 0];
r2_fix = [0 0];

r1 = [0 0];
r2 = [0 0];

angle = 90;

r_past = [0 0];
r1_past = [0 0];

% here is were the loop should be implemented

figure; grid; hold on

for i = 1:N-1

    % Implementation of offline minimum jerk methof
    A = [1; (10*(t/t_f)^3 -15*(t/t_f)^4 + 6*(t/t_f)^5)];
    X = [r_i' (r_f-r_i)'];
    
    r_off_list(i,:) = (X*A)';
    
    % The commented code will be used if online jerk has not to be used
    
    
    T = t_f - t; % time left
    tau = T/2;  % we are applying it every step, so the difference is dt
    
    % perturbation addition 
    T_mat = [T^3 T^4 T^5; 3*T^2 4*T^3 5*T^4; 6*T 12*T^2 20*T^3]; %3x3
    B = [(r_f-r-r1*T-(r2/2)*T^2); (-r1-r2*T); -r2]; % 3x2   I have to feed current state values
    A = T_mat\B; %3x2
      
    a0 = r;
    a1 = r1;
    a2 = r2/2;
    a3 = A(1,:);
    a4 = A(2,:);
    a5 = A(3,:);
    
    a = [a0' a1' a2' a3' a4' a5']; %2x6

    taus = [1; tau; tau^2; tau^3; tau^4; tau^5]; % 6x1
    
    taun = 0:dt:t_f-t;  
    tausn = [ones(1, size(taun,2)); taun; taun.^2; taun.^3; taun.^4; taun.^5]; % 6xN
    
    r_on_lists(i,:, 1:size(tausn,2)) = a*tausn;
    
    % update position
    r_on = a*taus;
    r_on = r_on';
       
    r =  r_off_list(i,:) + [r_on(1) 0];  % adding the x component predicted by online min jerk
    %pert = 0.001* rotateMat(r1_past, angle);
    r = r + [0.0000*norm(r1_past)  0];  % proportional to vel at previous time step
    %r = r + [pert(1) 0];  % proportional to vel at previous time step
    
    r1 = (r - r_past)/dt;
    r2 = (r1 - r1_past)/dt;
    
    % draw the trajectory
    %scatter(r(1), r(2), '.'); ylim([-5,15]); xlim([-10,10]);
    scatter(r_on(1), r_on(2), '.'); ylim([-5,15]); xlim([-10,10]);
    % draw the start and target points
    scatter(r_i(1), r_i(2), 'filled', 'k');
    scatter(r_f(1), r_f(2), 'filled', 'r');     

    r_past = r;
    r1_past = r1;
    r_list(i,:) = r;
    r_on_list(i,:) = r_on;
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




