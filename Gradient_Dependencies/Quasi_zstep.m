function [ gamma_step ] = Quasi_zstep(B_Quasi, G_n, )
J = diag([ones(xlength,1);-ones(ylength,1)]);
x1 = z_xi(:,1);
x2= z_xi(:,2);

y1= z_yj(:,1);
y2= z_yj(:,2);

gamma_step = - eta*B_Updated()*(G(z_xi,z_yj, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,  sigma, xlength, ylength));
end

