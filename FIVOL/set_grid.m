function [g] = set_grid(g)
%g is a structure containing grid information

for j = 2:g.nr-1
    jm = j-1;
    jp = j+1;
    for k = 2:g.nth-1
        km = k-1;
        kp = k+1;
        xa = 0.25*(g.x(j,k) + g.x(jm,k) + g.x(jm,km) + g.x(j,km));
        ya = 0.25*(g.y(j,k) + g.y(jm,k) + g.y(jm,km) + g.y(j,km));
        xb = 0.25*(g.x(j,k) + g.x(j,km) + g.x(jp,km) + g.x(jp,k));
        yb = 0.25*(g.y(j,k) + g.y(j,km) + g.y(jp,km) + g.y(jp,k));
        xc = 0.25*(g.x(j,k) + g.x(jp,k) + g.x(jp,kp) + g.x(j,kp));
        yc = 0.25*(g.y(j,k) + g.y(jp,k) + g.y(jp,kp) + g.y(j,kp));
        xd = 0.25*(g.x(j,k) + g.x(j,kp) + g.x(jm,kp) + g.x(jm,k));
        yd = 0.25*(g.y(j,k) + g.y(j,kp) + g.y(jm,kp) + g.y(jm,k));
        
        % Side AB
        
        dxa = xb - xa;
        dya = yb - ya;
        dxk = g.x(j,k) - g.x(j,km);
        dyk = g.y(j,k) - g.y(j,km);
        sab = abs(dxa*dyk - dxk*dya);
        g.Qab(j,k) = (dxa*dxa + dya*dya)/sab;
        g.Pab(j,k) = (dxa*dxk + dya*dyk)/sab;
        
        % Side BC
        
        dxb = xc - xb;
        dyb = yc - yb;
        dxj = g.x(j,k) - g.x(jp,k);
        dyj = g.y(j,k) - g.y(jp,k);
        sbc = abs(dyj*dxb - dxj*dyb);
        g.Qbc(j,k) = (dxb*dxb + dyb*dyb)/sbc;
        g.Pbc(j,k) = (dxb*dxj + dyb*dyj)/sbc;
        
        % Side CD
        
        dxc = xd - xc;
        dyc = yd - yc;
        dxk = g.x(j,k) - g.x(j,kp);
        dyk = g.y(j,k) - g.y(j,kp);
        scd = abs(dxc*dyk - dyc*dxk);
        g.Qcd(j,k) = (dxc*dxc + dyc*dyc)/scd;
        g.Pcd(j,k) = (dxc*dxk + dyc*dyk)/scd;
        
        % Side DA
        
        dxd = xa - xd;
        dyd = ya - yd;
        dxj = g.x(j,k) - g.x(jm,k);
        dyj = g.y(j,k) - g.y(jm,k);
        sda = abs(dxj*dyd - dyj*dxd);
        g.Qda(j,k) = (dxd*dxd + dyd*dyd)/sda;
        g.Pda(j,k) = (dxd*dxj + dyd*dyj)/sda;
    end
end

end

