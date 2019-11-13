function [ T ] = phi_grad(z, a1,a2,a3,a4,a5,a6,a7,a8,sigma)
x1=  z(:,1);
x2= z(:,2);



T0 = a1 + (1+a3)*x1 +a4 *x2 + a6 * (x1 -a7)*(-1/(sigma^2)).*exp(-0.5*(1/(sigma^2))*((x1-a7).^2 + (x2 -a8).^2));

T1 = a2 + (1+a5).*x2 +a4 .*x1 + a6 .* (x2 -a8).*(-1/(sigma^2)).*exp(-0.5*(1/(sigma^2))*((x1-a7).^2 + (x2 -a8).^2));

T = [T0 T1];




end
