function [ P , sig ] = optimizerOne(Meta, GeographicInfo, Wind, floodProbThreads, iter)
%
% Provide an optimize plan with given data
%
% input:
%   Meta.periodToBePlanned              how many periods do we plan on
%   Meta.maxIter                        what's the max of recursion
%   Meta.probIncrement 
%   Meta.kf             100
%   Meta.uf                     = 80;
%   Meta.carlength              = 80/3600*3; 
%
%   GeographicInfo.nodes = nodes
%   GeographicInfo.isMediator
%   GeographicInfo.ETable
%   GeographicInfo.ECap
%   GeographicInfo.ELength
%   GeographicInfo.ELanes
%   GeographicInfo.sizeV
%   GeographicInfo.population
%
% output:
%   P is the planner (for 
%

disp('Start planning iteration ');disp(iter);

if Meta.maxIter < iter || floodProbThreads >= 1
    disp('Exceed the maximum iterations allowed.');
    P = zeros(Meta.periodToBePlanned*GeographicInfo.sizeV,GeographicInfo.sizeV);
    sig = -1;
    return;
end;

%%
% Given the wind inofrmation, get the probablity of flooding
%
sizeV = GeographicInfo.sizeV;
probFlood = prob_flood( GeographicInfo.nodes , Wind);    % this is the probablity 
isSource = sign(max(0, probFlood - floodProbThreads));  % 0 if not flooded, 1 if flooded

%%
% solve for the first questions
% A*X<=b
A = [];
b = [];

% Aeq*X = beq
Aeq = [];
beq = [];

% lb <= X <= ub
% in this case  X >= 0
lb = zeros(size(GeographicInfo.ETable));
ub = [];

% use non linear constraint for it
% ceq   = @(X) [(sum((X.* GeographicInfo.ETable - (X.* GeographicInfo.ETable)'),2) .* isSource - GeographicInfo.population .* isSource);
%                (sum((X.* GeographicInfo.ETable - (X.* GeographicInfo.ETable)'),2) .* GeographicInfo.isMediator)];        
nonlcon  = @(X) deal([], [(sum((X.* GeographicInfo.ETable - (X.* GeographicInfo.ETable)'),2) .* isSource - GeographicInfo.population .* isSource);
               (sum((X.* GeographicInfo.ETable - (X.* GeographicInfo.ETable)'),2) .* GeographicInfo.isMediator)]);

% Create the cost functions
ELength = GeographicInfo.ELength;
ELanes  = GeographicInfo.ELanes;

% time    = @(r, c, fe) Meta.kf * (1+Meta.carlength*fe) * (Meta.carlength*fe + ELanes(r,c)*ELength(r,c)) / (Meta.kf * (ELanes(r,c)*Meta.uf*(1+Meta.carlength*fe) - Meta.uf * fe));

% generate array represents row and column indexes
c = 100;
cost = @(X) sum(sum(Meta.kf*(c*X.*X.*GeographicInfo.ETable + c*(ELanes.*ELength + 1) .* X.*GeographicInfo.ETable + ELanes.*ELength) ./ (1+ Meta.uf * (Meta.kf*c*ELanes - 1) .* (X + 1).*GeographicInfo.ETable + Meta.kf*Meta.uf*ELanes)));

% create the starting point
x0 = zeros(sizeV);

% generate options
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunEvals',60000);
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxFunEvals',600000);
[M,fval,exitflag,output] = fmincon(cost,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

%%
% clean up
M = M .* GeographicInfo.ETable;
fval

% check the flag
if exitflag > 1 % success!
    % find one of the path to check whether it exceeds the time limit (=4D)    
    repSourceIndex = find(isSource,1); 
    if isempty(repSourceIndex)
        % can return directory
        P = zeros(sizeV*Meta.periodToBePlanned,sizeV);
        return;
    end
    
    isSink = max(0, ones(sizeV,1) - isSource - GeographicInfo.isMediator);
    [ ~, s , hasCycle] = findPath( M, repSourceIndex, isSink);
    if hasCycle == 1
        disp('Cycles isnt expected in optimized solution.');
        P = zeros(sizeV*Meta.periodToBePlanned,sizeV);
        return;
    end
    disp('Time for theplan:');disp(s);
    
    if s < Meta.periodToBePlanned * 3600
        % this plan is good enough
        % spread out the plan according to the capacity of the roads
        P = zeros(sizeV*Meta.periodToBePlanned, sizeV);
        startIndex = Meta.periodToBePlanned - ceil(s);
        for t = 1 : Meta.periodToBePlanned
            if t >= startIndex % doesn't move unless need to
                P((t-1)*sizeV:t*sizeV, 1:sizeV) = min(M, GeographicInfo.ECap);
                M =  M - P((t-1)*sizeV:t*sizeV, 1:sizeV);
            else         
                P((t-1)*sizeV:t*sizeV, 1:sizeV) = zeros(sizeV); 
            end;
        end;
    else % then should recurse with higher threashold
        [P, sig] = optimizerOne(Meta, GeographicInfo, Wind, floodProbThreads+Meta.probIncrement, iter+1);
        if sig < 0
            P = M; % current plan could be the best
            sig = 1;
        end;
    end;
else
    % something is wrong that we cannot solve with the solver1 loop it up
    disp(M);disp(fval);disp(exitflag);disp(output);
    disp('Something is wrong with the optimizer');
    P = M;
    sig = exitflag;
    return
end;

end

