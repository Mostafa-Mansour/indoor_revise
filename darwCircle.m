
% problem number 4

function drawCircle(x,y,r)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);
axis equal
end

d

cdf_=@(x)pi/(sqrt(1-x^2)*x+asin(x));