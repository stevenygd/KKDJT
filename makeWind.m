%
% This secript will randomly generate input data
%
clear;
p = input('New input path:','s');
if ~exist('path','dir')
    mkdir(p);
end

%%
% Generate a wind
%
gen_seed = 0.01;
lat = input('lat:');
long = input('long:');
category = input('Input your category:');

center   = [(lat+rand()*gen_seed) (long+rand()*gen_seed)];

scale_to_cp = [
 % cate     base_cp     cp_range
    1       980         10000;
    2       965         (979-965);
    3       945         (964-945);
    4       920         (944-920);
    5       902         (920-902)
];

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
