%
% load the data
%   The simulator data includes:
%   1. V = nodes information        : |V| x 4   matrix
%       V_{i,1}             is identity, 
%       (V_{i,2},V_{i,3})   is the geographic point 
%       V_{i,4}             is the population
%
%   2. EC = edge capacity           : |V| x |V| matrix
%       EC_{i,j}            capacity of the edge from V_i to V_j
%
%   3. W  = Wind information        : 1 x       matrix
%       (W_{1}, W_{2})      is the geographic point for landfall
%       W_{3}               category
%

V = loadfile();

