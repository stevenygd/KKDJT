path = input('Path:','s');
P  = importdata(strcat(pwd, '/', path, '/plan.dat'));
N = size(P,2)
T = floor(size(P,1)/size(P,2))

[X, Y]= meshgrid(1:N, 1:N);

for t = 1:N
    D = P((1+(t-1)*N):(t*N),:);
    SP = sign(D(:,:));
    figure;
    title(strcat('Plan at Day'));
    spy(SP);
    for i = 1:N        
        hold on
        scatter(X(:,i), Y(:,i), (sqrt(Z(:,i))+10*log(Z(:,i))),'filled');
    end
    title(['Transportation Plan Day ',num2str(t)]);
    xlabel('Out Nodes ID')
    ylabel('In Nodes ID')
    grid on
    hold off
end
% 
% for t = 1:1
%     D = P((1+(t-1)*N):(t*N),:);
%     SP = sign(D(:,:));
%     Z = arrayfun(@(x,y) D(x,y), X, Y); 
%     figure; contour(Y, X, Z);
% end