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
P  = importdata(strcat(pwd, '/', path, '/prediction.dat'));

%
% Given the wind inofrmation, get the probablity of flooding
%
V               = N.data;
NodeTable       = N.textdata;
sizeV           = size(V, 1);
initPopulation  = V(:,3);

% extract edge parameters
ETable      = E(1:sizeV,:);
ECap        = E((sizeV+1):(sizeV*2),:);
ECostMove   = E((2*sizeV+1):(sizeV*3),:);

% use the wind information to calculate the prob flood


% solve for the first questions


time = @(edge, fe) (edge(1)+ fe^2)/(fe+edge(2));% TODO dummy!

