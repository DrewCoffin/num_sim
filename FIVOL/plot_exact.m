%plot exact solution to sin theta/r

n = 10

th = [0:pi/2/n:pi/2];
rxy = 2.0;
r0 = 0.5;
r = [r0:(rxy-r0)/n:rxy];
[TH,R] = meshgrid(th,r);

phi = sin(TH)./R;
plot(phi)
%h = polar([0 2*pi], [0 1]);
%delete(h)


[X,Y] = pol2cart(TH,R);

contour(X,Y,phi)

surf(X,Y,phi)
xlabel('x')
ylabel('y')
zlabel('\phi')
axis square