% 2 (more of this code was obtained from Prof. Smy)
clear

Boxes{1}.X = [0.8 1.2]*1e-7;
Boxes{1}.Y = [0.6 1.0]*1e-7;
Boxes{1}.BC = 0.0;

Boxes{2}.X = [0.8 1.2]*1e-7;
Boxes{2}.Y = [0.0 0.4]*1e-7;
Boxes{2}.BC = 0.0;

xlim = 200e-9;
ylim = 100e-9;
nx = 200;
ny = 100;
Vapp = 1.5;
Acond = 1;
Bcond = 0.0001;
q = 1.602e-19;  %Coul

[ Curr, Vmap, Ex, Ey, eFlowx, eFlowy  ] = ...
    Poisson(xlim,ylim,nx,ny,Acond,Bcond,[Vapp 0],Boxes);

figure(6)
subplot(2,1,1),H = surface(Vmap');
title('Voltage Map with Bottleneck')
set(H, 'linestyle', 'none');
view(45,45)
voltages=colorbar;
title(voltages,'Volts')

subplot(2,1,2)
quiver(eFlowx',eFlowy');
title('Electric Field from Potential - Vector Plot')
axis([0 nx 0 ny]);

%maxE = max(max(abs(eFlowx)))

Fx = q*eFlowx;
Fy = q*eFlowy;