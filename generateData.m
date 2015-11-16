%
% This secript will randomly generate input data
%
clear;
p = input('New input path:','s');
if ~exist('path','dir')
    mkdir(p);
end

%%
% save meta data
%
meta = [20505 0.5];
save(strcat(pwd, '/', p, '/meta.dat'), 'meta', '-ascii');

%%
% Generate a wind
%
gen_seed = 0.1;
center   = [(30.5+rand()*gen_seed) (89+rand()*gen_seed)];

scale_to_cp = [
 % cate     base_cp     cp_range
    1       74          (96-74);
    2       96          (111-96);
    3       111         (131-111);
    4       131         (115-131);
    5       955         (180-155)
];

category = input('Input your category:');
cp = scale_to_cp(category,2) + rand()*scale_to_cp(category,3);

wind = [category center(1) center(2) cp]
save(strcat(pwd, '/', p, '/wind.dat'), 'wind', '-ascii');

%%
% now make some prediction
%
prediction_error = [
 % time     cat_error   lat_error   long_error  cp_error
    1       0.0001      0.01        0.01        10      ;
    2       0.0001      0.01        0.01        3       ;
    3       0           0.005       0.005       1       ;
    4       0           0.001       0.001       0.5     ;
];

predictions = prediction_error(:,2:5) .* rand(4);
predictions = (1 + predictions) .* [wind; wind; wind; wind];
predictions(:,1) = floor(predictions(:,1));
predictions(2,:) = 0;
predictions
save(strcat(pwd, '/', p, '/prediction.dat'), 'predictions', '-ascii');

%%
% Generate some nodes
%
N = input('How many nodes do you want:');
node_range = [% approximate by mississipi map
    %   lat_base    long_base       population_base
        30          89              1000;
    %   lat_error   long_error      population_error
        0.8         2.0             242729
];
Nodes = [(1:N)' (rand(N,1)*node_range(2,1)+node_range(1,1)) (rand(N,1)*node_range(2,2)+node_range(1,2)) (rand(N,1)*node_range(2,3)+node_range(1,3)) sign(max(0,rand(N,1)-0.5))];
save(strcat(pwd, '/', p, '/nodes.dat'), 'Nodes', '-ascii');

%%
% Generate a graph
%
import Graph.*

capacity_cap = 81;
% generate edge capacity
E           = floor(rand(N)*sqrt(capacity_cap)); % -5 to reduce connectivity
E           = E .* sign(max(0,rand(N) - 0.9));
ECap        = E' * E;       % make it symmatrical since this is an undirected graph
ETable      = sign(ECap);

move_cap = 30;
E           = floor(rand(N)*sqrt(capacity_cap)); % -5 to reduce connectivity
ECostToMove = (E'*E).*ETable; % make it symmatrical since this is an undirected graph

length = 60;
ELength     = rand(N)*length + 30;

lanes = 4;
ELanes      = rand(N)*lanes + 2;

Edges = [ETable; ECap; ECostToMove; ELength; ELanes];
save(strcat(pwd, '/', p, '/edge.dat'), 'Edges', '-ascii');

figure
plot(graph(ETable, Nodes(:,1)));
figure
scatter(Nodes(:,2), Nodes(:,3), Nodes(:,4)./1000, 'filled');
