function [ pFlood ] = prob_flood( nodes )
%
% input:
%   nodes [sizeV x 4]        all the information about each nodes
%
% return :
%   pFlood sizeV x 1        probablity of flooding at each node
%
pFlood = 0.5 * magic(size(nodes,1));
pFlood = pFlood(:,1);
end

