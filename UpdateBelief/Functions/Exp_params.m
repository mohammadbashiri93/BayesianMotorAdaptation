function [ params ] = Exp_params( varargin )
% Exp_params returns the parameters for a specific experiment
% if no arguement is given the default experiment parameters are returned


% we initialize a struct array. each struct contains the parameters of a
% specific experiment:

% Example: params(1).name contains the name of the 1st experiment setup
%          params(2).rotMat containes the rotation function of the 
%          experiment 2

% I have initiated the function accepting variable number of arguements,
% maybe for later is useful. but for now we will limit the number of
% inputs.
errmsg = strcat('\nToo many arguements. The function only gets one input.',...
                '\n\nExample: \n \t exp = Exp_params(2)',...
                '\n\nThis saves the parameters of experiemnt 2 in exp');
assert (length(varargin) < 2, sprintf(errmsg));    
    

% Experiment 1 parameters
params(1).type    = 'CW';
params(1).compF   = @compForce1;


% Experiment 2 parameters
params(2).type    = 'CCW';
params(2).compF   = @compForce2;


ind = varargin{1};
if ~nargin
    % if no arguement was given return experiment 1's parameters
    params = params(1);   
else
    params = params(ind);
end

end

%% function to compute force given the velocity for each dynamic

% compute force applied for CW force field
function [force] = compForce1(v)
    k = 10;  
    RMat = [0 -1;  % angle = 90;
            1  0];
    force = k * v * RMat;
%     RMat = [cosd(angle) -sind(angle);
%             sind(angle) cosd(angle)]; % rotation matrix 
end

% compute force applied for CCW force field
function [force] = compForce2(v)
    k = 1;  
    RMat = [0  1;  % angle = -90;
            1  0];
    force = k * v * RMat;
%     RMat = [cosd(angle) -sind(angle);
%             sind(angle) cosd(angle)]; % rotation matrix 
end
