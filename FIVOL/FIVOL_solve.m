function [g,F] = FIVOL_solve(g, lambda, nmax, eps)
%g contains grid info
%lambda is the relaxation parameter
%nmax is the maximum number of iterations
%eps is the tolerance for convergence

surf(g.x,g.y,g.phi);
F = getframe;
rms = 1000.0;
drms = 1000.0;
eps = 0.001;
nit =0.0;

while (rms > eps)
%while (rms > eps) & (nit < nmax)
    %for n = 1:nmax
    sum = 0.0;
   for k = 2:g.nth-1
       km = k-1;
       kp = k+1;
       for j = 2:g.nr-1
           jm = j-1;
           jp = j+1;
           phd = 0.25*(g.Pcd(j,k) - g.Pda(j,k))*g.phi(jm,km) + ...
               (g.Qcd(j,k) + 0.25*(g.Pbc(j,k) - g.Pda(j,k)))*g.phi(j,kp) + ...
               0.25*(g.Pbc(j,k) - g.Pcd(j,k))*g.phi(jp,kp) + ...
               (g.Qda(j,k) + 0.25*(g.Pcd(j,k) - g.Pab(j,k)))*g.phi(jm,k) + ...
               (g.Qbc(j,k) + 0.25*(g.Pab(j,k) - g.Pcd(j,k)))*g.phi(jp,k) + ...
               0.25*(g.Pda(j,k) - g.Pab(j,k))*g.phi(jm,km) + ...
               (g.Qab(j,k) + 0.25*(g.Pda(j,k) - g.Pbc(j,k)))*g.phi(j,km) + ...
               0.25*(g.Pab(j,k) - g.Pbc(j,k))*g.phi(jp,km);
           phd = phd/(g.Qab(j,k) + g.Qbc(j,k) + g.Qcd(j,k) + g.Qda(j,k));
           dif = phd - g.phi(j,k);
           sum = sum + dif*dif;
           g.phi(j,k) = g.phi(j,k) + lambda*dif;                  
           %g.phi = boundary2(g.phi,g.nth,g.nr);
       end
   end
   surf(g.x,g.y,g.phi);
   F = [F,getframe];
   rms1 = rms;
   rms = sqrt(sum/((g.nr-2)*(g.nth-2)))
%   drms = abs(rms1-rms);
   nit = nit + 1
end

sprintf('Convergence achieved after %i steps',nit)

end

