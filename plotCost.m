path = input('path:','s');
cat = input('Category:');
OUT_COST  = importdata(strcat(pwd, '/', path, '/out_costs.dat'));
OUT_POPU  = importdata(strcat(pwd, '/', path, '/out_popInfo.dat'));
T = size(OUT_POPU,2);

figure
hold on
for t = 1 : T
     plot(sqrt(OUT_POPU(:,t)), ...
     'LineWidth',2,...
     'MarkerSize',10);    
 end
 title('Population Changes for Plan for Category 5');
 ylabel('Population (sqrt(people))');
 xlabel('Time Step (Day)');
 hold off

cmove   = OUT_COST(1,:);
figure;
plot(cmove);
title(['Moving Cost of Plan Category 5',num2str(cat)]);
xlabel('Time Step (Day)');
ylabel('Estimated Cost of Moving');

cecon   = OUT_COST(2,:);
figure;
plot(cecon);
title(['Economics Cost of Plan Category 5',num2str(cat)]);
xlabel('Time Step (Day)');
ylabel('Economics Opportunity Costs');

figure;
cdeath  = OUT_COST(3:size(OUT_COST,1),:) * 10^10;
plot(cdeath);
title(['Causulty Cost of Plan Category ',num2str(cat)]);
xlabel('Time Step (Day)');
ylabel('Cost caused by the Death');

