function [ r_out, v_out, a_out ] = updateMinJerk( r, r_actual, v, v_actual, a, F_forcefield, F_adapt, dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global m;

%     F_adapt = -0.8*F_forcefield;
    F_jerk  = m.*a;
    F_error = 25.*(r - r_actual)+ 5.*(v - v_actual);
    F_hand  = F_jerk + F_error + F_adapt;
    F_total = F_hand + F_forcefield;
    
    a_T = F_total./m;
    
    v_out = v_actual + a_T*dt;  % 1st equation of motion
    r_out = r_actual + v_actual*dt + a_T*dt^2;
    a_out = (v_out - v_actual)/dt;

end

