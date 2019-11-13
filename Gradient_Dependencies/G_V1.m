function [ Grad_final ] = G( z_xi,z_yj, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12, sigma, xlength, ylength)


f = phi_grad(z_xi, a1,a2,a3,a4,a5,a6,a7,a8 ,sigma);
w = grad_y_g(f, b1,b2,b3,b4,b5,b6,b7,b8,b9, b10,b11,b12, sigma);
p = grad_alpha_phi(z_xi, a1,a2,a3,a4,a5,a6,a7,a8 ,sigma);


p1= p(1:xlength,1:8);

p2= p((xlength+1):2*xlength,1:8);

w1 = w(1:xlength,1);

w2 = w(1:xlength,2);
% 
% size(w1)
% size(p1)
% size(p2)
% size(w2)
chain_rule_grad_alpha = p1'*w1 + p2'*w2;

S_a = (1/xlength)*(chain_rule_grad_alpha);
% size(S_a)


term_two = (exp(highlighter_g(z_yj,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12, sigma)))'*grad_beta_g(z_yj, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12, sigma);
chain_rule_grad_beta= (ones(1,xlength)*grad_beta_g(phi_grad(z_xi,a1,a2,a3,a4,a5,a6,a7,a8,sigma),b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,sigma)) - term_two;
S_b = (1/ylength)*chain_rule_grad_beta;
% size(S_b)
Grad_final = [S_a' S_b]';

end
