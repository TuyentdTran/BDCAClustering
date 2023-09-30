% Draw a circle in R2
function draw_circle (xcenter,ycenter,radius,color)
format long
theta = linspace(0,2*pi,200);
x = xcenter + radius * cos(theta);
y = ycenter + radius * sin(theta);
plot(x,y,'Color',color,'Linewidth',1.5);
end
