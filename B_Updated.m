function [B_updated] = B_Updated(G_nminus1, G_n, B_nminus1, J)
% This is the quasi implicit gradient descent algorithm
%   Detailed explanation goes here
alpha = (norm(J * G_n - B_nminus1*G_nminus1))^2 / ((G_nminus1)' * (J * G_n - B_nminus1*G_nminus1));

alpha_star = sign(alpha)* min(abs(alpha), norm(B_nminus1));
B_updated = B_nminus1 + alpha_star * ((J * G_n - B_nminus1*G_nminus1) *(J * G_n - B_nminus1*G_nminus1)' )/ (norm((J * G_n - B_nminus1*G_nminus1)))^2;

end

