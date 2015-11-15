function [ pFlood ] = prob_flood( nodes , windInfo)
%
% input:
%   nodes       [sizeV x 4]     all the information about each nodes
%   windInfo    [1x4]           basic wind information
%
% return :
%   pFlood sizeV x 1            probablity of flooding at each node
%

% pFlood = 0.5 * magic(size(nodes,1));
% pFlood = pFlood(:,1) ./ ((windInfo(1,2) - nodes(:,1)).^ 2 + (windInfo(1,3) - nodes(:,2)).^2);

locations = nodes(:,1:2);

% CP is the centrol buffer
cp          = windInfo(1,4);
deltaP      = 1013 - cp;
if deltaP > 48
    % if greater storm: 
    mu_rp       = 406.2 * deltaP^(-0.711);
    sigma_rp    = 187.7 * deltaP^(-0.711);

    mu_vf       = 6.6;
    sigma_vf    = 2.8;

    TL = betapdf(locations, 10.229, 11.747);

elseif deltaP > 31
    % if normal strom
    mu_rp       = 79.58 * deltaP^(-0.33);
    sigma_rp    = 36.78 * deltaP^(-0.33);

    mu_vf       = 5.5;
    sigma_vf    = 2.5;
    
    TL = normpdf(locations,-9.9,58.7);
    TL = max(0, TL + 90) - 90;      % truncate everything that is < -90
    TL = min(0, TL - 90) + 90;      % truncate everything that is > 90
else
    % no storm at all
    pFlood(1:size(locations,1)) = 0;
    return
end;

% alpha 5 is a normal distribution
CP(1:size(locations,1),1) = cp;
RP = lognpdf((1013 - CP),mu_rp,sigma_rp);  % magic number, fark field pressure
VF = lognpdf(TL,mu_vf,sigma_vf);
XD = normpdf(locations, 0,1);

pFlood = CP .* RP .* TL .* VF .* XD;

end

