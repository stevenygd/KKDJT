function [ d_rate ] = death_rate( nodes, time )
%
% input:
%   nodes   [sizeV x 4]     all the nodes information
%   time    float           the time after the landfall
% return:
%   d_rate  [sizeV x 1]     the rate of calsulty

% TODO: dummy!
d_rate = ones(size(nodes,1),1) / (exp(time) + 10);

end

