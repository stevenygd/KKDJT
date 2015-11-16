function [ output_args ] = dist_flood( windInfo )
%
% Create the distributions for floodings
%
% input:
%   nodes       [sizeV x 4]     all the information about each nodes
%   windInfo    [1x4]           basic wind information
%
% return :
%   pFlood sizeV x 1            probablity of flooding at each node
%

% CP is the central pressure
cp          = windInfo(1,4);
deltaP      = 1013 - cp;

if deltaP > 48
    % if greater storm: 
    mu_rp       = 406.2 * deltaP^(-0.711);
    sigma_rp    = 187.7 * deltaP^(-0.711);

    mu_vf       = 6.6;
    sigma_vf    = 2.8;
    
    TL = betarnd(10.229, 11.747);

elseif deltaP > 31
    % if normal strom
    mu_rp       = 79.58 * deltaP^(-0.33);
    sigma_rp    = 36.78 * deltaP^(-0.33);

    mu_vf       = 5.5;
    sigma_vf    = 2.5;
    
    TL = normrnd(locations,-9.9,58.7);
    TL = max(0, TL + 90) - 90;      % truncate everything that is < -90
    TL = min(0, TL - 90) + 90;      % truncate everything that is > 90

end;

% alpha 5 is a normal distribution
RP = lognpdf((1013 - cp),mu_rp,sigma_rp);  % magic number, fark field pressure
VF = lognpdf(TL,mu_vf,sigma_vf);
XD = normpdf(locations, 0,1);

pFlood = (cp * RP) .* TL .* VF .* XD; % magic number


end

