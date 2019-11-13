function [ L ] = L( z_x, z_y,xlength,ylength, a1,a2,a3,a4,a5,a6,a7,a8, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,sigma)

Precomposed_T = phi_grad(z_x, a1,a2,a3,a4,a5,a6,a7,a8 ,sigma);
Composed_T = highlighter_g(Precomposed_T,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12, sigma);

second_term = exp(highlighter_g(z_y, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12, sigma));

L = (1/xlength)*sum(Composed_T,'all') - (1/ylength) * sum(second_term, 'all');

end

