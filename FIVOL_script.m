%-------------------------------------------------------------------------
% FIVOL Fletcher 5.2.3
% Peter Delamere
% 4/3/15
%-------------------------------------------------------------------------

nth = 21;  %Number of points in the theta direction
nr = 21;   %Number of points in the radia direction

rmin = 0.1;  %Minimum r
rmax = 1.0;  %Maximum r

lambda = 1.5;   %Relaxation parameter
nmax = 1000.0;  %Maximum number of iteractions
eps = 0.0001;    %Tolerance for convergence for SOR

% initialize the grid

r = [rmin:(rmax-rmin)/(nr-1):rmax];
th = [0:pi/2/(nth-1):pi/2];
[x,y,TH,R] = get_grid(th,r);

phi = zeros(size(x));

% exact solution
figure(1);
phix = sin(TH)./R;
surf(x,y,phix);

% set boundary condtions

%phi = boundary2(phi,nth,nr);
phi = boundary(phi,phix,nth,nr);
surf(x,y,phi);

% grid related parameters
g = struct('nth', nth, 'nr', nr, 'x', x, 'y', y,...
   'Qab', zeros(size(x)), 'Pab', zeros(size(x)), 'Qbc', zeros(size(x)),...
   'Pbc', zeros(size(x)), 'Qcd', zeros(size(x)), 'Pcd', zeros(size(x)),...
   'Qda', zeros(size(x)), 'Pda', zeros(size(x)),...
   'phi', phi);

g = set_grid(g);

% Solve
figure(2);
[g,F] = FIVOL_solve(g, lambda, nmax, eps);

%movie(F,1)
figure(3)
surf(g.x,g.y,abs(g.phi-phix));

totdif2 = 0.0;
for j = 1:g.nr
    for k = 1:g.nth
        totdif2 = totdif2 + (g.phi(j,k) - phix(j,k)).^2;
    end
end
rms = sqrt((totdif2)./(g.nr.*g.nth))

%-------------------------------------------------------------------------
