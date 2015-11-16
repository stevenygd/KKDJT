rand = input('Is random (1/0)');
if rand == 1
    category = 4;
    if category > 2
        baseCp = 48;
    else
        baseCp = 31;
    end
    wind = [category (rand()*10) (rand()*10) (baseCp+rand()*17)];
    [wind(2) wind(3)]
    [X, Y]= meshgrid((wind(2)+(-1:0.05:1)'), (wind(3)+(-1:0.05:1)'));
    Z = arrayfun(@(X,Y) prob_flood([X Y], wind), X, Y); 
    figure; contourf(X,Y,Z);
else
    path = input('Input path:','s');
    wind  = importdata(strcat(pwd, '/', path, '/wind.dat'));
    [X, Y]= meshgrid((wind(2)+(-1:0.05:1)'), (wind(3)+(-1:0.05:1)'));
    Z = arrayfun(@(X,Y) prob_flood([X Y], wind), X, Y); 
    figure; contour(Y, X, Z);
    
    Nodes = importdata(strcat(pwd, '/', path, '/nodes.dat'));
    Edges = importdata(strcat(pwd, '/', path, '/edge.dat'));
    hold on
    scatter(Nodes(:,3), Nodes(:,2), Nodes(:,4)./1000,'filled');
    hold off
    figure; plot(graph(Edges(1:length(Nodes),:), Nodes(:,1)));

end;


