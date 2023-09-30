function []=plot_figurev2(C,R,CB,RB,optX,I)
hold on
[k,~]=size(optX);
qk=size(CB,1);
m=size(C,1);

color=[0.8500 0.3250 0.0980];
for i=1:m
    draw_circle(C(i,1),C(i,2),R(i),color);
end

color=[0 0.4470 0.7410];
for i = 1:k
    temp=projball2(optX(i,:),C(I==i,:),R(I==i));
    for j = 1:size(temp,1)
        plot([optX(i,1),temp(j,1)],[optX(i,2),temp(j,2)],'Color',color,'Linewidth',1.5)
    end
end

plot(optX(:,1),optX(:,2),'bp','MarkerFaceColor','b','MarkerSize',15)

for i=1:qk
    draw_circle(CB(i,1),CB(i,2),RB(i),'black');
end
end