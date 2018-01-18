function [ sspace ] = retrieve_memory( varargin )
% MEMORY acts as the implicit memory


errmsg = strcat('\nToo many arguements. The function only gets one input.',...
                '\n\nExample: \n \t exp = memory(2)',...
                '\n\nThis saves the parameters of experiemnt 2 in exp');
assert (length(varargin) < 2, sprintf(errmsg));

% this should preferrably be an input to the function!
global Vx;
global Vy;

Vsize = numel(Vx); % TODO: must come from Globals!!! DONE!

memory(1).name    = 'null';
memory(1).Fmus    = zeros(2, Vsize, Vsize);
memory(1).Fsigmas = ones(2, Vsize, Vsize)*10;

memory(2).name    = 'CW';
memory(2).Fmus    = zeros(2, Vsize, Vsize);
memory(2).Fsigmas = ones(2, Vsize, Vsize)*10;

memory(3).name    = 'CCW';
memory(3).Fmus    = zeros(2, Vsize, Vsize);
memory(3).Fsigmas = ones(2, Vsize, Vsize)*10;

% ...

ind = varargin{1};
if ~nargin
    % if no arguement was given return null field memory
    sspace = memory(1);   
else
    sspace = memory(ind);
end

end
