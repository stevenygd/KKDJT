%
% load the data
%   The simulator data includes:
%   0.  M = meta data                : 1 x 2
%       M(1,1) GDP per capita
%       M(1,2) threashhold for the flooding probablity
%
%   1.  V = nodes information        : |V|  x 4   matrix
%       V_{i,1}             is identity, 
%       (V_{i,2},V_{i,3})   is the geographic point 
%       V_{i,4}             is the population
%
%   2.  E = edge parameters          : |V|  x |V| matrix       
%       First   |V|x|V| matrix       : adjencency matrix
%       Second  |V|x|V| matrix       : capactiy of of each edge
%                                       (0 if not connected)
%       Third  |V|x|V| matrix       : cost of moving alone each edge
%                                       (Infinity if not connected)
%
%   3.  W  = Wind information        : 1    x 4   matrix
%       W(1)               is the time period of the landfall
%       (W(2), W(3))       is the geographic point for landfall
%       W(4)               category
%       W(5)               central pressure
%
%   4.  P  = Plan information        : n|V| x |V| matrix
%       each |V|x|V| matrix denotes a plan day for a period
%       plan for the first n period;
%       if n > [numPeriod] then ignore all the next
%
%   Output:
%       simulatedPeriods = p        how many time period is simualted
%       cmove            1xp        the cost to move at each time period
%       cecon            1xp        the opp cost at each time period
%       cdeath           px|V|      the expected calsulty at each period
%                                   for each node           
%       popInfo          (p+1)x|V|  the remaining population for each time
%                                   stamp and each node            
%
%

%
% Load Data
%
import Graph.*
pwd

path = input('Input path:','s');
numPeriod = input('Number of simulation period:');

M  = importdata(strcat(pwd, '/', path, '/meta.dat'));
N  = importdata(strcat(pwd, '/', path, '/nodes.dat'));
E  = importdata(strcat(pwd, '/', path, '/edge.dat'));
W  = importdata(strcat(pwd, '/', path, '/wind.dat'));
P  = importdata(strcat(pwd, '/', path, '/plan.dat'));

%
% Simulate
%
V = N(:,2:4);
% NodeTable = N.textdata;
sizeV           = size(V, 1);
initPopulation  = V(:,3);

plannedPeriods  = size(P,1) / sizeV;

landfallTime    = W(1,1);

gpaPerCapita    = M(1,1);
floodCap        = M(1,2);

% extract edge parameters
ETable      = E(1:sizeV,:);
ECap        = E((sizeV+1):(sizeV*2),:);
ECostMove   = E((2*sizeV+1):(sizeV*3),:);

% step for each day
simulatedPeriods = max(max(plannedPeriods, numPeriod), landfallTime + 2);
cmove   = zeros(1, simulatedPeriods);  % [0xn] mat. contains costs
                                       % to move at each simulated period
                                                    
cecon   = zeros(1, simulatedPeriods);  % [0xn] mat. contains opportunity cost
                                       % at each simulated period
                                       
% [n x |V|] mat. contains the casulty at each period 
% should expect this to be zero before the landfall
cdeath  = zeros(simulatedPeriods, sizeV);  

% [1+n, V] mat contains population remaining at each county at each time period
popInfo = zeros(simulatedPeriods + 1, sizeV);
popInfo(1,:) = initPopulation';
                                                    
for t = 1 : simulatedPeriods
    if t > plannedPeriods
        p = zeros(sizeV, sizeV) .* ETable;
    else 
        p = P((1+sizeV*(t-1)):sizeV*t,:) .* ETable;
    end
    
    % calculate the cost to move at each period
    pExceed = arrayfun(@max,zeros(sizeV),p - ECap);
    costExceedMat = costMatrixOfExceed(pExceed);
    cmove(1,t) = sum(sum(p .* ECostMove + costExceedMat));
    
    % calculate the economic opportunity cost at each period
    totalMove  = sum(sum(p));
    op_cost    = @(diff) max(0,diff) ^ 2;
    cecon(1,t) = op_cost(landfallTime - t) * totalMove * gpaPerCapita;
    
    % update the population
    vOut = sum(p,2)';    % out population for each vertices  (sum alone rows)
    vIn  = sum(p);       % in population for each verices    (sum alone cols)
    popInfo(t+1,:) = max(0, popInfo(t,:) - vOut + vIn);
    
    % death cost
    if t >= landfallTime
       % TODO: floating point error! consider!
       deathForEachNode = sign(max(0, prob_flood(V, W) - floodCap)) .* death_rate(V, t - landfallTime);
       cdeath(t,:) = deathForEachNode' .* popInfo(t,:);       

       % disp('Death incurred by the wind:');disp(t);disp(deathForEachNode);
       
       popInfo(t+1,:) = max(0, popInfo(t,:) .* (1- deathForEachNode'));
    end;
    
    % disp(popInfo(t+1,:))
end;

% cmove
cecon
% popInfo
% cdeath

save(strcat(pwd, '/', path, '/out_costs.dat'), 'cmove',   'cecon', 'cdeath', '-ascii','-tabs');
save(strcat(pwd, '/', path, '/out_popInfo.dat'), 'popInfo', '-ascii','-tabs'); 


