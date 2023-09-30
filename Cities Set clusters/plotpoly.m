function plotpoly(Z,col)
%% plot ploygons
kp=length(Z);
hold on
for j=1:kp
    z=Z{j};
    %fill([z(:,1);z(1,1)],[z(:,2);z(1,2)],col)
    plot([z(:,1);z(1,1)],[z(:,2);z(1,2)],'Color',col,'LineWidth',1.5)
end
hold off
end