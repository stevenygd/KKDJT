%
% load the data
%   The simulator data includes:
%   0.  M = meta data                : 1 x 2
%       M(1,1) GDP per capita
%       M(1,2) threashhold for the flooding probablity
%
%   1.  V = nodes information        : |V|  x 4   matrix
%       V(i,1)             identity, 
%       (V(i,2),V(i,3))    the geographic point 
%       V(i,4)             the population
%       V(i,5)             mediator nodes
%
%   2.  E = edge parameters          : |V|  x |V| matrix       
%       First   |V|x|V| matrix       : adjencency matrix
%       Second  |V|x|V| matrix       : capactiy of of each edge
%                                       (0 if not connected)
%       Third  |V|x|V| matrix       : cost of moving alone each edge
%                                       (Infinity if not connected)
%
%   3.  W  = Wind information        : 1    x 4   matrix
%       W_{1}               is the time period of the landfall
%       (W_{2}, W_{3})      is the geographic point for landfall
%       W_{4}               category
%
%   4.  P  = Prediction              : 4 x ?    matrix
%       P(i,:) is the prediction for the ith timestamp
%       P(2,:) is always zero
%
%   Output:
%   Plan : 
%
clear;
path = input('Input path:','s');
floodProbThreads = input('Threashold for the probability:');
totalPeriodsToBePlanned = input('How many periods in total should we planned:');

M  = importdata(strcat(pwd, '/', path, '/meta.dat'));
N  = importdata(strcat(pwd, '/', path, '/nodes.dat'));
E  = importdata(strcat(pwd, '/', path, '/edge.dat'));
W  = importdata(strcat(pwd, '/', path, '/prediction.dat'));

%%
% Set up the geographic data
V               = N(:,2:5);
% NodeTable       = N.textdata;
initPopulation  = V(:,3);
isMediator      = V(:,4);
sizeV           = size(V, 1);

HasEdge     = E(1:sizeV,:);
ECap        = E((sizeV+1):(sizeV*2),:);
ECostMove   = E((2*sizeV+1):(sizeV*3),:);
ELength     = E((3*sizeV+1):(sizeV*4),:);
ELanes      = E((4*sizeV+1):(sizeV*5),:);

GeographicInfo.nodes        = V;
GeographicInfo.isMediator   = isMediator;
GeographicInfo.ETable       = HasEdge;
GeographicInfo.ECap         = ECap;
GeographicInfo.ECostMove    = ECostMove;
GeographicInfo.V            = V;
% GeographicInfo.NodeTable    = NodeTable;
GeographicInfo.sizeV        = sizeV;
GeographicInfo.population   = initPopulation;
GeographicInfo.ELength      = ELength;
GeographicInfo.ELanes       = ELanes;

%%
% Set up the prediction data
%
Meta.maxIter                = 10;
Meta.probIncrement          = 0.03;
Meta.kf                     = 100;
Meta.uf                     = 80;
Meta.carlength              = 80/3600*3; 

P = zeros(totalPeriodsToBePlanned * sizeV, sizeV);      % store the actual plan
TP = zeros(totalPeriodsToBePlanned * sizeV, sizeV);     % store the tmp plan
for t = 1:totalPeriodsToBePlanned
    if size(W,1) < t % has no more information, copy the rest
        numLeft = totalPeriodsToBePlanned - t + 1;
        P((sizeV*(t-1)+1):sizeV*(numLeft+t-1),:) = TP(sizeV:sizeV*(numLeft + 1), :);     
        break;
    end
    
    windInfo = W(1,:);
    % do a new planning if has the new information isn't void 
    % i.e. the wind information is not [0 0 0 ... 0]
    if find(windInfo) > 0
        Meta.periodToBePlanned      = 4 - t - 1;

        % set up for optimizer one
        [ TP ] = optimizerOne(Meta, GeographicInfo, windInfo, 0.5, 1);
    end
    
    % update the geographic information
    CurrPlan = TP(1:sizeV,:)
    P((sizeV*(t-1)+1):sizeV*t,:) = CurrPlan;
    GeographicInfo.population = GeographicInfo.population + sum(CurrPlan,2) - sum(CurrPlan,2);
    GeographicInfo.population'
end

P
save(strcat(pwd, '/', path, '/plan.dat'), 'P', '-ascii');










