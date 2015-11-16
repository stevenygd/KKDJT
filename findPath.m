function [ path, s , hasCycle] = findPath( E, i, isSink )
% find the path from the ith node to one of the nodes in Sink
% 
% input:
%   E           |V|x|V|     edge graphs.
%   isSource    |V|x1       bit map, whether it is source.
%   isSink      |V|x1       bit map, whether it is sink.
%

if sum(sum(isSink)) == 0 || i < 1 || i > length(E) 
    disp('Error:No souce or sinks.');
    return
end

path = zeros(1, length(E)*2-1); % 1x(len of longest path I expected)
count = 1;
path(count) = i;
s = 0;
curr = i;
E = E.* (ones(length(E)) - diag(ones(1,length(E)))); % clean up, put the diagnal to zero
hasCycle = 0;
while isSink(curr,1) ~= 1 % if curr isn't a sink
  newCurr = find(E(curr,:),1);
  
  s     = s + E(curr, newCurr);  
  curr  = newCurr;  
  count = count + 1;
  path(1,count) = curr;
  
  % detach cycle
  if sum(ismember(path, newCurr)) > 1      
      hasCycle = 1;
      disp('Detect cycle in the graph.');
      return;
  end;
end

end

