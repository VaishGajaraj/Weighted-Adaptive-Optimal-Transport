function [grad_alpha_phi] = grad_alpha_phi(z, a1,a2,a3,a4,a5,a6,a7,a8 ,sigma)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

x=  z(:,1);
y= z(:,2);
exp_arg = exp(-0.5*(1/(sigma^2))*((x-a7).^2 + (y -a8).^2));


alpha_1_5= (x - a7).*(-1/(2*sigma^2)).*exp_arg;
alpha_2_5 = (y - a8).*(-1/(2*sigma^2)).*exp_arg;
alpha_1_6= -a6*(-1/(sigma^2)).*exp_arg + a6.*(x-a7).*(1/2*sigma^4).*(-2*x-2*a7) .*exp_arg;
alpha_2_6= a6.*(y-a8).*(1/2*sigma^4).*(-2*x-2*a7) .*exp_arg;
alpha_1_7= a6.*(x-a6).*(1/2*sigma^4).*(-2*y-2*a8) .*exp_arg;
alpha_2_7 = -a6.*(-1/(sigma^2)).*exp_arg + a6.*(y-a8).*(1/2*sigma^4).*(-2*y-2*a8) .* exp_arg;

first_row= [ones(size(x)) zeros(size(x)) x y zeros(size(x)) alpha_1_5 alpha_1_6 alpha_1_7];
second_row = [zeros(size(x)) ones(size(x)) zeros(size(x)) x y alpha_2_5 alpha_2_6 alpha_2_7];

grad_alpha_phi = [first_row; second_row];

end

