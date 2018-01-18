function [ r, v, a ] = compMinJerk( startPos, targetPos, period, t)
%COMPMINJERK Summary of this function goes here
%   Detailed explanation goes here
    
    tau    = t/period;
    r = startPos + (startPos-targetPos) .* ...
        (15*tau^4 - 6*tau^5 - 10*tau^3);
    v = (startPos-targetPos) .* (60*tau^3/period - ...
        30*tau^4/period - 30*tau^2/period);
    a = (startPos-targetPos) .* (180*tau^2/period^2 - ...
        120*tau^3/period^2 - 60*tau/period^2);

end

